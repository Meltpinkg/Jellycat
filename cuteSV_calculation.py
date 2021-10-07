import numpy as np
import math

'''
    node -> list(Record)
    return represent_idx, (center)start, end
'''
def cal_center(node):
    #return len(node) / 2
    r_start = 0
    r_end = 0
    start_list = []
    end_list = []
    for record in node:
        r_start += record.start
        r_end += record.end
        start_list.append(record.start)
        end_list.append(record.end)
    r_start = r_start / len(node)
    r_end = r_end / len(node)
    r_id = 0
    start_dif = 0x3f3f3f3f
    for id in range(len(node)):
        if abs(node[id].start - r_start) < start_dif:
            start_dif = abs(node[id].start - r_start)
            r_id = id
    return r_id, cal_ci(start_list), cal_ci(end_list)

def cal_ci(input_list):
    pos = int(1.96 * np.std(input_list) / len(input_list) ** 0.5)
    return "-%d,%d"%(pos, pos)

'''
    semi cluster for ins
    node -> dict{id -> list(Record)}
    candidates -> [list(Record)] -> ([clu1, clu2, clu3, ...])
    clu -> list(Record)
    UPDATE 2021.08.29
    check whether merged by check_two_sv_INS
    merge each two sample by km
'''
def cal_can_INS(node, debug):
    if debug == 1:
        print('start cal can:')
        for id in node:
            print(str(id) + ': ')
            for r in node[id]:
                print(r.to_string())
    candidates = list()  # clusters
    can_len = 0
    can_idx = 0
    for id in node:
        if len(node[id]) > can_len:
            can_idx = id
            can_len = len(node[id])

    for record in node[can_idx]:
        candidates.append([{can_idx}, record.start, record.end, [record], 1])
    for id in node:
        if id == can_idx:
            continue
        if len(candidates) == 1:
            if len(node[id]) > 1:
                print('ERROR when len(candidates) == 1') #####
            if debug:
                print(id)
                print(candidates)
                print(node[id][0].to_string())
            if check_two_sv_INS(candidates[0][1], candidates[0][2], node[id][0].start, node[id][0].end):
                candidates[0][0].add(id)
                candidates[0][1] = (candidates[0][4] * candidates[0][1] + node[id][0].start) / (candidates[0][4] + 1)
                candidates[0][2] = (candidates[0][4] * candidates[0][2] + node[id][0].end) / (candidates[0][4] + 1)
                candidates[0][3].append(node[id][0])
                candidates[0][4] += 1
            else:
                candidates.append([{id}, node[id][0].start, node[id][0].end, [node[id][0]], 1])
        elif len(candidates) < 100:
            if debug == 1:
                print(candidates)
                for record in node[id]:
                    print(record.to_string())
            dis = []
            n = len(candidates)
            for record in node[id]:
                dis.append([])
                for candidate in candidates:
                    if check_two_sv_INS(candidate[1], candidate[2], record.start, record.end):
                        distance = -int( math.log(abs(candidate[1] - record.start) + 1) * 100 + math.log(abs(candidate[2] - record.end) + 1) * 1000 )
                        #distance = -math.log(abs(candidate[2] - record.end) + 1)
                        '''if min(candidate[2], record.end) == 0:
                            bias = abs(max(candidate[2], record.end) * 10) + 1
                            print(bias)
                        else:
                            bias = abs(max(candidate[2], record.end) / min(candidate[2], record.end) * 10) + 1'''
                    else:
                        distance = -0x3f3f3f3f
                    dis[-1].append(distance)
            for i in range(len(node[id]), n, 1):
                dis.append([-0x3f3f3f3f for j in range(n)])
            if debug:
                print(dis)
            linker = km(n, dis)
            if debug:
                for i in range(n):
                    print('%d %d: %d'%(linker[i], i, dis[linker[i]][i]))
            # merge
            for i in range(n):
                if dis[linker[i]][i] == -0x3f3f3f3f:
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                elif linker[i] == -1:
                    print('ERROR linker[i] == -1') #####
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                else:
                    record = node[id][linker[i]]
                    candidates[i][0].add(id)
                    candidates[i][1] = (candidates[i][4] * candidates[i][1] + record.start) / (candidates[i][4] + 1)
                    candidates[i][2] = (candidates[i][4] * candidates[i][2] + record.end) / (candidates[i][4] + 1)
                    candidates[i][3].append(record)
                    candidates[i][4] += 1
            if debug == 1:
                print(candidates)
            '''
            #boys = [(363, 509)]
            #girls = [(887, 51), (399, 509), (1012, 75)]
            '''
        else:
            #print('find len(candidates) >= 100: %s, go greedy strategy'%(candidates[0][3][0].type))
            st_bias = [0 for i in range(len(candidates))]
            ed_bias = [0 for i in range(len(candidates))]
            merge_flag = False
            for i in range(len(candidates)):
                if id not in candidates[i][0]:
                    if check_two_sv_INS(candidates[i][1], candidates[i][2], record.start, record.end):
                        merge_flag = True
                    st_bias[i] = abs(candidates[i][1] - record.start)
                    ed_bias[i] = abs(candidates[i][2] - record.end)
                else:
                    st_bias[i] = -1
                    ed_bias[i] = -1
            if debug == 1:
                print(merge_flag)
            if merge_flag is False:
                candidates.append([{id}, record.start, record.end, [record], 1])
                continue
            st_mean = cal_mean(st_bias)
            ed_mean = cal_mean(ed_bias)
            if ed_mean == 0 or st_mean == 0:
                ratio = 1
            else:
                ratio = st_mean / ed_mean
            if debug == 1:
                print(st_bias)
                print(ed_bias)
                print(ratio)
            ratio = 1 ###
            min_bias = 0x3f3f3f3f
            min_bias_idx = 0
            '''
            for i in range(len(st_bias)):
                if debug == 1:
                    print(st_bias[i] / ratio + ed_bias[i])
                if st_bias[i] == -1:
                    continue
                if st_bias[i] / ratio + ed_bias[i] < min_bias:
                    min_bias = st_bias[i] / ratio + ed_bias[i]
                    min_bias_idx = i
            '''
            for i in range(len(st_bias)):
                if debug == 1:
                    print(st_bias[i] / ratio + ed_bias[i])
                if st_bias[i] == -1:
                    continue
                st_bias[i] += 1
                ed_bias[i] += 1
                if math.log(st_bias[i]) + math.log(ed_bias[i]) < min_bias:
                    min_bias = math.log(st_bias[i]) + math.log(ed_bias[i])
                    min_bias_idx = i
            candidates[min_bias_idx][0].add(id)
            candidates[min_bias_idx][1] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][1] + record.start) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][2] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][2] + record.end) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][3].append(record)
            candidates[min_bias_idx][4] += 1
    ans_can = list()
    for candidate in candidates:
        ans_can.append(candidate[3])
    # re-merge candidates
    '''
    ans_can = list()
    len_can = len(candidates)
    vis_id = set()
    for i in range(len_can):
        if i in vis_id:
            continue
        vis_id.add(i)
        id_temp = candidates[i][0]
        st_temp = candidates[i][1]
        ed_temp = candidates[i][2]
        record_temp = candidates[i][3]
        for j in range(i, len_can, 1):
            if j not in vis_id:
                if len(id_temp & candidates[j][0]) == 0 and check_two_sv_INS(candidates[i][1], candidates[i][2], candidates[j][1], candidates[j][2]):
                    id_temp = id_temp | candidates[j][0]
                    st_temp = (st_temp + candidates[j][1]) / 2
                    ed_temp = (ed_temp + candidates[j][2]) / 2
                    record_temp = record_temp + candidates[j][3]
                    vis_id.add(j)
        ans_can.append(record_temp)
    '''
    if debug:
        print('finish cal can')
        for i in range(len(ans_can)):
            print(str(i + 1) + ':')
            for r in ans_can[i]:
                print(r.to_string())
    return ans_can

def cal_can_DEL(node, debug):
    out_km = False
    if debug == 1:
        print('start cal can:')
        for id in node:
            print(str(id) + ': ')
            for r in node[id]:
                print(r.to_string())
    candidates = list()  # clusters
    can_len = 0
    can_idx = 0
    for id in node:
        if len(node[id]) > can_len:
            can_idx = id
            can_len = len(node[id])

    for record in node[can_idx]:
        candidates.append([{can_idx}, record.start, record.end, [record], 1])
    for id in node:
        if id == can_idx:
            continue
        if len(candidates) == 1:
            if len(node[id]) > 1:
                print('ERROR when len(candidates) == 1') #####
            if debug:
                print(id)
                print(candidates)
                print(node[id][0].to_string())
            if check_two_sv_DEL(candidates[0][1], candidates[0][2], node[id][0].start, node[id][0].end):
                candidates[0][0].add(id)
                candidates[0][1] = (candidates[0][4] * candidates[0][1] + node[id][0].start) / (candidates[0][4] + 1)
                candidates[0][2] = (candidates[0][4] * candidates[0][2] + node[id][0].end) / (candidates[0][4] + 1)
                candidates[0][3].append(node[id][0])
                candidates[0][4] += 1
            else:
                candidates.append([{id}, node[id][0].start, node[id][0].end, [node[id][0]], 1])
        elif len(candidates) < 100:
            if debug == 1:
                print(candidates)
                for record in node[id]:
                    print(record.to_string())
            dis = []
            n = len(candidates)
            for record in node[id]:
                dis.append([])
                for candidate in candidates:
                    if check_two_sv_DEL(candidate[1], candidate[2], record.start, record.end):
                        bias = abs((max(candidate[2], record.end) + 1) / (min(candidate[2], record.end) + 1) * 10)
                        distance = -int( math.log(abs(candidate[1] - record.start) + 1) * 100 + math.log(bias) * 1000 )
                    else:
                        distance = -0x3f3f3f3f
                    dis[-1].append(distance)
            for i in range(len(node[id]), n, 1):
                dis.append([-0x3f3f3f3f for j in range(n)])
            if debug:
                print(dis)
            linker = km(n, dis)
            if debug:
                for i in range(n):
                    print('%d %d: %d'%(linker[i], i, dis[linker[i]][i]))
            # merge
            for i in range(n):
                if dis[linker[i]][i] == -0x3f3f3f3f:
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                elif linker[i] == -1:
                    print('ERROR linker[i] == -1') #####
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                else:
                    record = node[id][linker[i]]
                    candidates[i][0].add(id)
                    candidates[i][1] = (candidates[i][4] * candidates[i][1] + record.start) / (candidates[i][4] + 1)
                    candidates[i][2] = (candidates[i][4] * candidates[i][2] + record.end) / (candidates[i][4] + 1)
                    candidates[i][3].append(record)
                    candidates[i][4] += 1
            if debug == 1:
                print(candidates)
        else:
            out_km = True
            st_bias = [0 for i in range(len(candidates))]
            ed_bias = [0 for i in range(len(candidates))]
            merge_flag = False
            for i in range(len(candidates)):
                if id not in candidates[i][0]:
                    if check_two_sv_DEL(candidates[i][1], candidates[i][2], record.start, record.end):
                        merge_flag = True
                    st_bias[i] = abs(candidates[i][1] - record.start)
                    ed_bias[i] = abs(candidates[i][2] - record.end)
                else:
                    st_bias[i] = -1
                    ed_bias[i] = -1
            if debug == 1:
                print(merge_flag)
            if merge_flag is False:
                candidates.append([{id}, record.start, record.end, [record], 1])
                continue
            st_mean = cal_mean(st_bias)
            ed_mean = cal_mean(ed_bias)
            if ed_mean == 0 or st_mean == 0:
                ratio = 1
            else:
                ratio = st_mean / ed_mean
            if debug == 1:
                print(st_bias)
                print(ed_bias)
                print(ratio)
            ratio = 1 ###
            min_bias = 0x3f3f3f3f
            min_bias_idx = 0
            for i in range(len(st_bias)):
                if debug == 1:
                    print(st_bias[i] / ratio + ed_bias[i])
                if st_bias[i] == -1:
                    continue
                st_bias[i] += 1
                ed_bias[i] += 1
                if math.log(st_bias[i]) + math.log(ed_bias[i]) < min_bias:
                    min_bias = math.log(st_bias[i]) + math.log(ed_bias[i])
                    min_bias_idx = i
            candidates[min_bias_idx][0].add(id)
            candidates[min_bias_idx][1] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][1] + record.start) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][2] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][2] + record.end) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][3].append(record)
            candidates[min_bias_idx][4] += 1
    ans_can = list()
    for candidate in candidates:
        ans_can.append(candidate[3])
    if debug:
        print('finish cal can')
        for i in range(len(ans_can)):
            print(str(i + 1) + ':')
            for r in ans_can[i]:
                print(r.to_string())
    if out_km:
        print('ERROR len(candidates) >= 100: %s'%(candidates[0][3][0].type)) #####
        #for cc in candidates:
        #    print(cc[1],end='')
        #print()
    return ans_can

def cal_can_INV(node, debug):
    out_km = False
    if debug == 1:
        print('start cal can:')
        for id in node:
            print(str(id) + ': ')
            for r in node[id]:
                print(r.to_string())
    candidates = list()  # clusters
    can_len = 0
    can_idx = 0
    for id in node:
        if len(node[id]) > can_len:
            can_idx = id
            can_len = len(node[id])

    for record in node[can_idx]:
        candidates.append([{can_idx}, record.start, record.end, [record], 1])
    for id in node:
        if id == can_idx:
            continue
        if len(candidates) == 1:
            if len(node[id]) > 1:
                print('ERROR when len(candidates) == 1') #####
            if debug:
                print(id)
                print(candidates)
                print(node[id][0].to_string())
            if check_two_sv_INV(candidates[0][1], candidates[0][2], node[id][0].start, node[id][0].end):
                candidates[0][0].add(id)
                candidates[0][1] = (candidates[0][4] * candidates[0][1] + node[id][0].start) / (candidates[0][4] + 1)
                candidates[0][2] = (candidates[0][4] * candidates[0][2] + node[id][0].end) / (candidates[0][4] + 1)
                candidates[0][3].append(node[id][0])
                candidates[0][4] += 1
            else:
                candidates.append([{id}, node[id][0].start, node[id][0].end, [node[id][0]], 1])
        elif len(candidates) < 100:
            if debug == 1:
                print(candidates)
                for record in node[id]:
                    print(record.to_string())
            dis = []
            n = len(candidates)
            for record in node[id]:
                dis.append([])
                for candidate in candidates:
                    if check_two_sv_INV(candidate[1], candidate[2], record.start, record.end):
                        distance = -abs(candidate[1] - record.start) - abs(candidate[2] - record.end) - 1
                    else:
                        distance = -0x3f3f3f3f
                    dis[-1].append(distance)
            for i in range(len(node[id]), n, 1):
                dis.append([-0x3f3f3f3f for j in range(n)])
            if debug:
                print(dis)
            linker = km(n, dis)
            if debug:
                for i in range(n):
                    print('%d %d: %d'%(linker[i], i, dis[linker[i]][i]))
            # merge
            for i in range(n):
                if dis[linker[i]][i] == -0x3f3f3f3f:
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                elif linker[i] == -1:
                    print('ERROR linker[i] == -1') #####
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                else:
                    record = node[id][linker[i]]
                    candidates[i][0].add(id)
                    candidates[i][1] = (candidates[i][4] * candidates[i][1] + record.start) / (candidates[i][4] + 1)
                    candidates[i][2] = (candidates[i][4] * candidates[i][2] + record.end) / (candidates[i][4] + 1)
                    candidates[i][3].append(record)
                    candidates[i][4] += 1
            if debug == 1:
                print(candidates)
        else:
            out_km = True
            st_bias = [0 for i in range(len(candidates))]
            ed_bias = [0 for i in range(len(candidates))]
            merge_flag = False
            for i in range(len(candidates)):
                if id not in candidates[i][0]:
                    if check_two_sv_INV(candidates[i][1], candidates[i][2], record.start, record.end):
                        merge_flag = True
                    st_bias[i] = abs(candidates[i][1] - record.start)
                    ed_bias[i] = abs(candidates[i][2] - record.end)
                else:
                    st_bias[i] = -1
                    ed_bias[i] = -1
            if debug == 1:
                print(merge_flag)
            if merge_flag is False:
                candidates.append([{id}, record.start, record.end, [record], 1])
                continue
            st_mean = cal_mean(st_bias)
            ed_mean = cal_mean(ed_bias)
            if ed_mean == 0 or st_mean == 0:
                ratio = 1
            else:
                ratio = st_mean / ed_mean
            if debug == 1:
                print(st_bias)
                print(ed_bias)
                print(ratio)
            ratio = 1 ###
            min_bias = 0x3f3f3f3f
            min_bias_idx = 0
            for i in range(len(st_bias)):
                if debug == 1:
                    print(st_bias[i] / ratio + ed_bias[i])
                if st_bias[i] == -1:
                    continue
                st_bias[i] += 1
                ed_bias[i] += 1
                if math.log(st_bias[i]) + math.log(ed_bias[i]) < min_bias:
                    min_bias = math.log(st_bias[i]) + math.log(ed_bias[i])
                    min_bias_idx = i
            candidates[min_bias_idx][0].add(id)
            candidates[min_bias_idx][1] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][1] + record.start) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][2] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][2] + record.end) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][3].append(record)
            candidates[min_bias_idx][4] += 1
    ans_can = list()
    for candidate in candidates:
        ans_can.append(candidate[3])
    if debug:
        print('finish cal can')
        for i in range(len(ans_can)):
            print(str(i + 1) + ':')
            for r in ans_can[i]:
                print(r.to_string())
    if out_km:
        print('find len(candidates) >= 100: %s, go greedy strategy'%(candidates[0][3][0].type)) #####
        #for cc in candidates:
        #    print(cc[1],end='')
        #print()
    return ans_can

def cal_can_DUP(node, debug):
    if debug == 1:
        print('start cal can:')
        for id in node:
            print(str(id) + ': ')
            for r in node[id]:
                print(r.to_string())
    candidates = list()  # clusters
    can_len = 0
    can_idx = 0
    for id in node:
        if len(node[id]) > can_len:
            can_idx = id
            can_len = len(node[id])

    for record in node[can_idx]:
        candidates.append([{can_idx}, record.start, record.end, [record], 1])
    for id in node:
        if id == can_idx:
            continue
        if len(candidates) == 1:
            if len(node[id]) > 1:
                print('ERROR when len(candidates) == 1') #####
            if debug:
                print(id)
                print(candidates)
                print(node[id][0].to_string())
            if check_two_sv_DUP(candidates[0][1], candidates[0][2], node[id][0].start, node[id][0].end):
                candidates[0][0].add(id)
                candidates[0][1] = (candidates[0][4] * candidates[0][1] + node[id][0].start) / (candidates[0][4] + 1)
                candidates[0][2] = (candidates[0][4] * candidates[0][2] + node[id][0].end) / (candidates[0][4] + 1)
                candidates[0][3].append(node[id][0])
                candidates[0][4] += 1
            else:
                candidates.append([{id}, node[id][0].start, node[id][0].end, [node[id][0]], 1])
        elif len(candidates) < 100:
            if debug == 1:
                print(candidates)
                for record in node[id]:
                    print(record.to_string())
            dis = []
            n = len(candidates)
            for record in node[id]:
                dis.append([])
                for candidate in candidates:
                    if check_two_sv_DUP(candidate[1], candidate[2], record.start, record.end):
                        #bias = abs((max(candidate[2], record.end) + 1) / (min(candidate[2], record.end) + 1) * 10)
                        #distance = -int( math.log(abs(candidate[1] - record.start) + 1) * 100 + math.log(bias) * 1000 )
                        distance = -abs(candidate[1] - record.start) - abs(candidate[2] - record.end) - 1
                    else:
                        distance = -0x3f3f3f3f
                    dis[-1].append(distance)
            for i in range(len(node[id]), n, 1):
                dis.append([-0x3f3f3f3f for j in range(n)])
            if debug:
                print(dis)
            linker = km(n, dis)
            if debug:
                for i in range(n):
                    print('%d %d: %d'%(linker[i], i, dis[linker[i]][i]))
            # merge
            for i in range(n):
                if dis[linker[i]][i] == -0x3f3f3f3f:
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                elif linker[i] == -1:
                    print('ERROR linker[i] == -1') #####
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                else:
                    record = node[id][linker[i]]
                    candidates[i][0].add(id)
                    candidates[i][1] = (candidates[i][4] * candidates[i][1] + record.start) / (candidates[i][4] + 1)
                    candidates[i][2] = (candidates[i][4] * candidates[i][2] + record.end) / (candidates[i][4] + 1)
                    candidates[i][3].append(record)
                    candidates[i][4] += 1
            if debug == 1:
                print(candidates)
            '''
            #boys = [(363, 509)]
            #girls = [(887, 51), (399, 509), (1012, 75)]
            '''
        else:
            #print('ERROR len(candidates) >= 100: %s'%(candidates[0][3][0].type)) #####
            st_bias = [0 for i in range(len(candidates))]
            ed_bias = [0 for i in range(len(candidates))]
            merge_flag = False
            for i in range(len(candidates)):
                if id not in candidates[i][0]:
                    if check_two_sv_DUP(candidates[i][1], candidates[i][2], record.start, record.end):
                        merge_flag = True
                    st_bias[i] = abs(candidates[i][1] - record.start)
                    ed_bias[i] = abs(candidates[i][2] - record.end)
                else:
                    st_bias[i] = -1
                    ed_bias[i] = -1
            if debug == 1:
                print(merge_flag)
            if merge_flag is False:
                candidates.append([{id}, record.start, record.end, [record], 1])
                continue
            st_mean = cal_mean(st_bias)
            ed_mean = cal_mean(ed_bias)
            if ed_mean == 0 or st_mean == 0:
                ratio = 1
            else:
                ratio = st_mean / ed_mean
            if debug == 1:
                print(st_bias)
                print(ed_bias)
                print(ratio)
            ratio = 1 ###
            min_bias = 0x3f3f3f3f
            min_bias_idx = 0
            for i in range(len(st_bias)):
                if debug == 1:
                    print(st_bias[i] / ratio + ed_bias[i])
                if st_bias[i] == -1:
                    continue
                st_bias[i] += 1
                ed_bias[i] += 1
                if math.log(st_bias[i]) + math.log(ed_bias[i]) < min_bias:
                    min_bias = math.log(st_bias[i]) + math.log(ed_bias[i])
                    min_bias_idx = i
            candidates[min_bias_idx][0].add(id)
            candidates[min_bias_idx][1] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][1] + record.start) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][2] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][2] + record.end) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][3].append(record)
            candidates[min_bias_idx][4] += 1
    ans_can = list()
    for candidate in candidates:
        ans_can.append(candidate[3])
    if debug:
        print('finish cal can')
        for i in range(len(ans_can)):
            print(str(i + 1) + ':')
            for r in ans_can[i]:
                print(r.to_string())
    return ans_can

def cal_can_BND(node, debug):
    if debug == 1:
        print('start cal can:')
        for id in node:
            print(str(id) + ': ')
            for r in node[id]:
                print(r.to_string())
    candidates = list()  # clusters
    can_len = 0
    can_idx = 0
    for id in node:
        if len(node[id]) > can_len:
            can_idx = id
            can_len = len(node[id])

    for record in node[can_idx]:
        candidates.append([{can_idx}, record.start, record.end, [record], 1])
    for id in node:
        if id == can_idx:
            continue
        if len(candidates) == 1:
            if len(node[id]) > 1:
                print('ERROR when len(candidates) == 1') #####
            if debug:
                print(id)
                print(candidates)
                print(node[id][0].to_string())
            if check_two_sv_BND(candidates[0][1], candidates[0][2], node[id][0].start, node[id][0].end):
                candidates[0][0].add(id)
                candidates[0][1] = (candidates[0][4] * candidates[0][1] + node[id][0].start) / (candidates[0][4] + 1)
                candidates[0][2] = (candidates[0][4] * candidates[0][2] + node[id][0].end) / (candidates[0][4] + 1)
                candidates[0][3].append(node[id][0])
                candidates[0][4] += 1
            else:
                candidates.append([{id}, node[id][0].start, node[id][0].end, [node[id][0]], 1])
        elif len(candidates) < 100:
            if debug == 1:
                print(candidates)
                for record in node[id]:
                    print(record.to_string())
            dis = []
            n = len(candidates)
            for record in node[id]:
                dis.append([])
                for candidate in candidates:
                    if check_two_sv_BND(candidate[1], candidate[2], record.start, record.end):
                        bias = abs((max(candidate[2], record.end) + 1) / (min(candidate[2], record.end) + 1) * 10)
                        distance = -int( math.log(abs(candidate[1] - record.start) + 1) * 100 + math.log(bias) * 1000 )
                    else:
                        distance = -0x3f3f3f3f
                    dis[-1].append(distance)
            for i in range(len(node[id]), n, 1):
                dis.append([-0x3f3f3f3f for j in range(n)])
            if debug:
                print(dis)
            linker = km(n, dis)
            if debug:
                for i in range(n):
                    print('%d %d: %d'%(linker[i], i, dis[linker[i]][i]))
            # merge
            for i in range(n):
                if dis[linker[i]][i] == -0x3f3f3f3f:
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                elif linker[i] == -1:
                    print('ERROR linker[i] == -1') #####
                    if linker[i] < len(node[id]):
                        candidates.append([{id}, node[id][linker[i]].start, node[id][linker[i]].end, [node[id][linker[i]]], 1])
                else:
                    record = node[id][linker[i]]
                    candidates[i][0].add(id)
                    candidates[i][1] = (candidates[i][4] * candidates[i][1] + record.start) / (candidates[i][4] + 1)
                    candidates[i][2] = (candidates[i][4] * candidates[i][2] + record.end) / (candidates[i][4] + 1)
                    candidates[i][3].append(record)
                    candidates[i][4] += 1
            if debug == 1:
                print(candidates)
            '''
            #boys = [(363, 509)]
            #girls = [(887, 51), (399, 509), (1012, 75)]
            '''
        else:
            #print('ERROR len(candidates) >= 100: %s'%(candidates[0][3][0].type)) #####
            st_bias = [0 for i in range(len(candidates))]
            ed_bias = [0 for i in range(len(candidates))]
            merge_flag = False
            for i in range(len(candidates)):
                if id not in candidates[i][0]:
                    if check_two_sv_BND(candidates[i][1], candidates[i][2], record.start, record.end):
                        merge_flag = True
                    st_bias[i] = abs(candidates[i][1] - record.start)
                    ed_bias[i] = abs(candidates[i][2] - record.end)
                else:
                    st_bias[i] = -1
                    ed_bias[i] = -1
            if debug == 1:
                print(merge_flag)
            if merge_flag is False:
                candidates.append([{id}, record.start, record.end, [record], 1])
                continue
            st_mean = cal_mean(st_bias)
            ed_mean = cal_mean(ed_bias)
            if ed_mean == 0 or st_mean == 0:
                ratio = 1
            else:
                ratio = st_mean / ed_mean
            if debug == 1:
                print(st_bias)
                print(ed_bias)
                print(ratio)
            ratio = 1 ###
            min_bias = 0x3f3f3f3f
            min_bias_idx = 0
            for i in range(len(st_bias)):
                if debug == 1:
                    print(st_bias[i] / ratio + ed_bias[i])
                if st_bias[i] == -1:
                    continue
                st_bias[i] += 1
                ed_bias[i] += 1
                if math.log(st_bias[i]) + math.log(ed_bias[i]) < min_bias:
                    min_bias = math.log(st_bias[i]) + math.log(ed_bias[i])
                    min_bias_idx = i
            candidates[min_bias_idx][0].add(id)
            candidates[min_bias_idx][1] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][1] + record.start) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][2] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][2] + record.end) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][3].append(record)
            candidates[min_bias_idx][4] += 1
    ans_can = list()
    for candidate in candidates:
        ans_can.append(candidate[3])
    if debug:
        print('finish cal can')
        for i in range(len(ans_can)):
            print(str(i + 1) + ':')
            for r in ans_can[i]:
                print(r.to_string())
    return ans_can

def check_two_sv_INS(start1, end1, start2, end2):
    max_ins_dist = 1000
    max_ins_ratio = 0.3
    mean = max(end1, end2)
    if mean == 0:
        end1 = 1
        end2 = 1
    if abs(start1 - start2) < max_ins_dist:
        if mean < 100 and min(end1, end2) / max(end1, end2) > max_ins_ratio:
            return True
        if 100 <= mean < 500 and min(end1, end2) / max(end1, end2) > max_ins_ratio:
            return True
        if 500 <= mean < 1000 and min(end1, end2) / max(end1, end2) > max_ins_ratio:
            return True
        if mean >= 1000 and min(end1, end2) / max(end1, end2) > max_ins_ratio:
            return True
    elif abs(start1 - start2) < max_ins_dist * 1.5:
        if min(end1, end2) / max(end1, end2) > 0.9:
            return True
    return False

def check_two_sv_DEL(start1, end1, start2, end2):
    max_dist = 1000
    mean = max(end1, end2)
    if mean == 0:
        end1 = 1
        end2 = 1
    if abs(start1 - start2) < max_dist:
        if mean < 100 and min(end1, end2) / max(end1, end2) > 0.3:
            return True
        if 100 <= mean < 500 and min(end1, end2) / max(end1, end2) > 0.4:
            return True
        if 500 <= mean < 1000 and min(end1, end2) / max(end1, end2) > 0.2:
            return True
        if mean >= 1000 and min(end1, end2) / max(end1, end2) > 0.2:
            return True
    elif abs(start1 - start2) < max_dist * 1.5:
        if min(end1, end2) / max(end1, end2) > 0.9:
            return True
    return False

def check_two_sv_INV(start1, end1, start2, end2):
    max_dist = 1000
    mean = max(end1, end2)
    if mean == 0:
        end1 = 1
        end2 = 1
    #end1 = start1 + end1
    #end2 = start2 + end2
    if abs(start1 - start2) < max_dist and abs(end1 - end2) < 300:
        return True
    return False

def check_two_sv_DUP(start1, end1, start2, end2):
    max_dist = 1000
    mean = max(end1, end2)
    if mean == 0:
        end1 = 1
        end2 = 1
    if abs(start1 - start2) < max_dist:
        if min(end1, end2) / max(end1, end2) > 0.3:
            return True
    elif abs(start1 - start2) < max_dist * 1.5:
        if min(end1, end2) / max(end1, end2) > 0.9:
            return True
    return False

def check_two_sv_BND(start1, end1, start2, end2):
    max_dist = 800
    if abs(start1 - start2) < max_dist and abs(end1 - end2) < max_dist:
        return True
    return False

def cal_mean(bias):
    sum = 0
    num = 0
    for i in bias:
        if i != -1:
            sum += i
            num += 1
    return sum / num

def dfs(x, n, visx, visy, lx, ly, dis, linker, slack):
    visx[x] = 1
    for y in range(n):
        if visy[y]:
            continue
        tmp = lx[x] + ly[y] - dis[x][y]
        if tmp < 1e-6:
            visy[y] = 1
            if linker[y] == -1 or dfs(linker[y], n, visx, visy, lx, ly, dis, linker, slack):
                linker[y] = x
                return 1
        elif slack[y] > tmp:
            slack[y] = tmp
    return 0


def km(n, dis):
    linker = [-1 for i in range(n)]
    ly = [0 for i in range(n)]
    lx = [-0x3f3f3f3f for i in range(n)]
    for i in range(n):
        for j in range(n):
            if dis[i][j] > lx[i]:
                lx[i] = dis[i][j]
    #print(lx)
    for x in range(n):
        slack = [0x3f3f3f3f for i in range(n)]
        while True:
            visx = [0 for i in range(n)]
            visy = [0 for i in range(n)]
            if dfs(x, n, visx, visy, lx, ly, dis, linker, slack):
                break
            d = 0x3f3f3f3f
            for i in range(n):
                if (not visy[i]) and d > slack[i]:
                    d = slack[i]
            for i in range(n):
                if visx[i]:
                    lx[i] -= d
            for i in range(n):
                if visy[i]:
                    ly[i] += d
                else:
                    slack[i] -= d
    return linker