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

def pre_forest(node, svtype, debug):
    points = []
    idx = 0
    edges = []
    for id in node:
        for i in range(len(node[id])):
            record = node[id][i]
            points.append([id, i, record.start, record.end])
            idx += 1
    # coordinate x
    points = sorted(points, key = lambda x:x[2])
    for i in range(len(points)):
        for j in range(i, min(i + 10, len(points))):
            if i != j:
                if check_two_sv(svtype, points[i][2], points[i][3], points[j][2], points[j][3]) and points[i][0] != points[j][0]:
                    distance = -cal_distance(svtype, points[i][2], points[j][2], points[i][3], points[j][3])
                    #print(points[i])
                    #print(points[j])
                    #print(distance)
                    edges.append([points[i][0], points[i][1], points[j][0], points[j][1], distance])
    # coordinate y
    points = sorted(points, key = lambda x:x[3])
    #print(len(points))
    for i in range(len(points)):
        for j in range(i, min(i + 5, len(points))):
            if i != j:
                if check_two_sv(svtype, points[i][2], points[i][3], points[j][2], points[j][3]) and points[i][0] != points[j][0]:
                    distance = -cal_distance(svtype, points[i][2], points[j][2], points[i][3], points[j][3])
                    edges.append([points[i][0], points[i][1], points[j][0], points[j][1], distance])
    edges = sorted(edges, key = lambda x:x[4])
    new_edges = []
    uniq_edges = []
    i = 0
    #[[1, 1, 2, 1, 69], [1, 1, 2, 1, 69], [1, 0, 2, 0, 4122], [2, 0, 1, 0, 4122], [0, 0, 1, 0, 4924], [1, 0, 0, 0, 4924], [0, 0, 2, 0, 5350], [2, 0, 0, 0, 5350]]
    while i < len(edges):
        new_edges.append(edges[i])
        while i+1 < len(edges) and ((edges[i][0] == edges[i+1][2] and edges[i][1] == edges[i+1][3] and edges[i][2] == edges[i+1][0] and edges[i][3] == edges[i+1][1]) or \
            (edges[i][0] == edges[i+1][0] and edges[i][1] == edges[i+1][1] and edges[i][2] == edges[i+1][2] and edges[i][3] == edges[i+1][3])):
            i += 1
        i += 1
    preline = {}
    for edge in new_edges:
        if edge[0] != edge[2] and (edge[0], edge[1]) not in preline and (edge[2], edge[3]) not in preline:
            preline[(edge[0], edge[1])] = (edge[2], edge[3])
            preline[(edge[2], edge[3])] = (edge[0], edge[1])
        if len(preline) == min(len(new_edges) / 4, 5) * 2:
            break
    if debug:
        print('edges')
        print(edges)
        print('preline:')
        print(preline)
    return preline

def merge_forest(preline, candidates, records, id, svtype, debug):
    import time
    st = time.time()
    if debug == 1:
        print('start merging forest')
        print(candidates)
        for record in records:
            print(record.to_string())
    dis = []
    new_candidates = []
    remove_candidates = set()
    remove_records = set()
    for i in range(len(records)):
        if (id, i) in preline:
            break_flag1 = False
            for can_i in range(len(candidates)):
                if can_i in remove_candidates:
                    continue
                candidate = candidates[can_i]
                for pr in candidate[5]:
                    if pr == preline[(id, i)]:
                        candidate[0].add(id)
                        candidate[1] = (candidate[4] * candidate[1] + records[i].start) / (candidate[4] + 1)
                        candidate[2] = (candidate[4] * candidate[2] + records[i].end) / (candidate[4] + 1)
                        candidate[3].append(records[i])
                        candidate[4] += 1
                        candidate[5].append((id, i))
                        new_candidates.append(candidate)
                        break_flag1 = True
                        remove_candidates.add(can_i)
                        remove_records.add(i)
                        break
                if break_flag1:
                    break
    rest_candidates = []
    rest_records = []
    trans_records = []
    record_i = 0
    for i in range(len(candidates)):
        if i not in remove_candidates:
            rest_candidates.append(candidates[i])
    for i in range(len(records)):
        if i not in remove_records:
            rest_records.append(records[i])
            trans_records.append(i)
            record_i += 1
    if debug:
        print('new:')
        print(new_candidates)
        print(rest_candidates)
        print(rest_records)
    n = len(rest_candidates)
    #print('m1: %s '%(time.time() - st),end='')
    st = time.time()
    for record in rest_records:
        dis.append([])
        for candidate in rest_candidates:
            if check_two_sv(svtype, candidate[1], candidate[2], record.start, record.end):
                distance = cal_distance(svtype, candidate[1], record.start, candidate[2], record.end)
            else:
                distance = -0x3f3f3f3f
            dis[-1].append(distance)
    for i in range(len(rest_records), n, 1):
        dis.append([-0x3f3f3f3f for j in range(n)])
    if debug == 1:
        print(dis)
    linker = km(n, dis)
    if debug == 1:
        for i in range(n):
            print('%d %d: %d'%(linker[i], i, dis[linker[i]][i]))
    # merge
    for i in range(n):
        if dis[linker[i]][i] == -0x3f3f3f3f:
            if linker[i] < len(rest_records):
                rest_candidates.append([{id}, rest_records[linker[i]].start, rest_records[linker[i]].end, [rest_records[linker[i]]], 1, [(id, trans_records[linker[i]])]])
        elif linker[i] == -1:
            print('ERROR linker[i] == -1') #####
            if linker[i] < len(rest_records):
                rest_candidates.append([{id}, rest_records[linker[i]].start, rest_records[linker[i]].end, [rest_records[linker[i]]], 1, [(id, trans_records[linker[i]])]])
        else:
            record = rest_records[linker[i]]
            rest_candidates[i][0].add(id)
            rest_candidates[i][1] = (rest_candidates[i][4] * rest_candidates[i][1] + record.start) / (rest_candidates[i][4] + 1)
            rest_candidates[i][2] = (rest_candidates[i][4] * rest_candidates[i][2] + record.end) / (rest_candidates[i][4] + 1)
            rest_candidates[i][3].append(record)
            rest_candidates[i][4] += 1
            rest_candidates[i][5].append((id, trans_records[linker[i]]))
    if debug == 1:
        print('rest_candidates')
        print(rest_candidates)
    #print('m2: %s '%(time.time() - st),end='')
    return new_candidates + rest_candidates

'''
    semi cluster
    node -> dict{id -> list(Record)}
    candidates -> [list(Record)] -> ([clu1, clu2, clu3, ...])
    clu -> list(Record)
    UPDATE 2021.08.29
    check whether merged by check_two_sv
    merge each two sample by km
'''
def cal_can(node, debug, svtype):
    import time
    start_time = time.time()
    pretime = -1
    out_km = False
    gen_forest = True
    if debug == 1:
        print('==========start cal can:')
        for id in node:
            print(str(id) + ': ')
            for r in node[id]:
                print(r.to_string())
    candidates = list()  # clusters
    can_len = 0
    can_idx = 0
    id_length = dict()
    for id in node:
        id_length[id] = len(node[id])
        if len(node[id]) > can_len:
            can_idx = id
            can_len = len(node[id])
    id_length = sorted(id_length.items(), key = lambda x:x[1], reverse=True)
    for i in range(len(node[can_idx])):
        record = node[can_idx][i]
        candidates.append([{can_idx}, record.start, record.end, [record], 1, [(can_idx, i)]])
    #for id in node:
    for idd in id_length:
        id = idd[0]
        if id == can_idx:
            continue
        '''
        if len(candidates) == 1:
            if len(node[id]) > 1:
                print('ERROR when len(candidates) == 1') #####
            if debug:
                print(id)
                print(candidates)
                print(node[id][0].to_string())
            if check_two_sv(svtype, candidates[0][1], candidates[0][2], node[id][0].start, node[id][0].end):
                candidates[0][0].add(id)
                candidates[0][1] = (candidates[0][4] * candidates[0][1] + node[id][0].start) / (candidates[0][4] + 1)
                candidates[0][2] = (candidates[0][4] * candidates[0][2] + node[id][0].end) / (candidates[0][4] + 1)
                candidates[0][3].append(node[id][0])
                candidates[0][4] += 1
            else:
                candidates.append([{id}, node[id][0].start, node[id][0].end, [node[id][0]], 1])
        '''
        if len(candidates) < 100:
            '''
            if gen_forest:
                t1 = time.time()
                preline = pre_forest(node, svtype, debug)
                gen_forest = False
                pretime = time.time() - t1
                #print('pre_forest: %f  '%(time.time() - start_time),end='')
            candidates = merge_forest(preline, candidates, node[id], id, svtype, debug)
            '''
            if debug == 1:
                print(candidates)
                for record in node[id]:
                    print(record.to_string())
            dis = []
            n = len(candidates)
            for record in node[id]:
                dis.append([])
                for candidate in candidates:
                    if check_two_sv(svtype, candidate[1], candidate[2], record.start, record.end):
                        distance = cal_distance(svtype, candidate[1], record.start, candidate[2], record.end)
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
                    #print('ERROR linker[i] == -1') #####
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
            #'''
            '''
            #boys = [(363, 509)]
            #girls = [(887, 51), (399, 509), (1012, 75)]
            '''
        else:
            out_km = True
            for record in node[id]:
                st_bias = [0 for i in range(len(candidates))]
                ed_bias = [0 for i in range(len(candidates))]
                merge_flag = False
                for i in range(len(candidates)):
                    if id not in candidates[i][0]:
                        if check_two_sv(svtype, candidates[i][1], candidates[i][2], record.start, record.end):
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
                if len(id_temp & candidates[j][0]) == 0 and check_two_sv(svtype, candidates[i][1], candidates[i][2], candidates[j][1], candidates[j][2]):
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
    if out_km and debug:
        print('find len(candidates) >= 100: %s, go greedy strategy'%(candidates[0][3][0].type))
        for cc in candidates:
            print(cc[0], end='')
            for rr in cc[3]:
                print(', ' + rr.to_string(), end='')
            print(",%d,%d "%(int(cc[1]), int(cc[2])),end='')
        print()
    '''endtime = time.time()
    if pretime == -1:
        print((endtime-start_time)*100000)
    else:
        print('%f %f %f'%(pretime*100000, (endtime-start_time-pretime)*100000, (endtime-start_time)*100000))'''
    return ans_can

# distance = cal_distance(svtype, candidate[1], record.start, candidate[2], record.end)
def cal_distance(svtype, start1, start2, end1, end2):
    if svtype == 'DEL':
        bias = abs((max(end1, end2) + 1) / (min(end1, end2) + 1) * 10)
        return -int( math.log(abs(start1 - start2) + 1) * 100 + math.log(bias) * 1000 )
    if svtype == 'INV':
        return -abs(start1 - start2) - abs(end1 - end2) - 1
    if svtype == 'DUP':
        return -abs(start1 - start2) - abs(end1 - end2) - 1
    if svtype == 'BND':
        bias = abs((max(end1, end2) + 1) / (min(end1, end2) + 1) * 10)
        return -int( math.log(abs(start1 - start2) + 1) * 100 + math.log(bias) * 1000 )
    #if svtype == 'INS':
    else:
        return -int( math.log(abs(start1 - start2) + 1) * 100 + math.log(abs(end1 - end2) + 1) * 1000 )
    

def check_two_sv(svtype, start1, end1, start2, end2):
    if svtype == 'DEL':
        if check_two_sv_DEL(start1, end1, start2, end2):
            return True
    elif svtype == 'INV':
        if check_two_sv_INV(start1, end1, start2, end2):
            return True
    elif svtype == 'DUP':
        if check_two_sv_DUP(start1, end1, start2, end2):
            return True
    elif svtype == 'BND':
        if check_two_sv_BND(start1, end1, start2, end2):
            return True
    else:
        if check_two_sv_INS(start1, end1, start2, end2):
            return True
    return False

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

'''
if len(candidates) == 1:
    if len(node[id]) > 1:
        print('ERROR when len(candidates) == 1') #####
    if debug:
        print(id)
        print(candidates)
        print(node[id][0].to_string())
    if check_two_sv(svtype, candidates[0][1], candidates[0][2], node[id][0].start, node[id][0].end):
        candidates[0][0].add(id)
        candidates[0][1] = (candidates[0][4] * candidates[0][1] + node[id][0].start) / (candidates[0][4] + 1)
        candidates[0][2] = (candidates[0][4] * candidates[0][2] + node[id][0].end) / (candidates[0][4] + 1)
        candidates[0][3].append(node[id][0])
        candidates[0][4] += 1
    else:
        candidates.append([{id}, node[id][0].start, node[id][0].end, [node[id][0]], 1])
'''