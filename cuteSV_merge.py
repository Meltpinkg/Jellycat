#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||  :  |||// \
#                / _||||| -:- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  ___/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#         \  \ `_.   \_ __\ /__ _/   .-` /  /
#     =====`-.____`.___ \_____/___.-`___.-'=====
#                       `=---='
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#               佛祖保佑         永无BUG
#####################################################
from cuteSV_input import *
from cuteSV_calculation import *
from cuteSV_output import output_result, solve_annotation, generate_header
from multiprocessing import Pool
from heapq import *
import os
import sys
import time
import argparse

def parse_para(config_file):
    para = []
    with open(config_file) as f:
        for line in f:
            seq = line.strip().split(' ')
            try:
                para.append(int(seq[1]))
            except:
                para.append(float(seq[1]))
    return para
def parse_config_para(p, svtype):
    if svtype == 'DEL':
        return [p[1], p[7], p[8]]
    if svtype == 'INV':
        return [p[2], p[9], p[10]]
    if svtype == 'DUP':
        return [p[3], p[11], p[12]]
    if svtype == 'BND':
        return [p[4], p[13], p[14]]
    else: # INS and other SV type
        return [p[0], p[5], p[6]]

# vcf_files->dict{fileidx->list(Record)}  only one svtype
# return list[]
def first_sort(vcf_files, max_dist, filenum, x):
    heap = []
    ans = []
    start_down = 0
    start_up = 0
    ratio = 1.6
    record_pointer = [0 for i in range(filenum)]
    #print(vcf_files.keys())
    for fileidx in vcf_files:
        heappush(heap, vcf_files[fileidx][0])
    while len(heap) > 0:
        head = heappop(heap) # Record
        if len(ans) == 0 or (abs(head.start - start_down) > max_dist * ratio and abs(head.start - start_up) > max_dist * ratio):
            #ans.append([head])
            ans.append(dict())
            ans[-1][head.source] = [head]
            start_down = head.start
            start_up = head.start
        else:
            #ans[-1].append(head)
            if head.source in ans[-1]:
                ans[-1][head.source].append(head)
            else:
                ans[-1][head.source] = [head]
            start_down = min(start_down, head.start)
            start_up = max(start_up, head.start)
        record_pointer[head.source] += 1
        if record_pointer[head.source] < len(vcf_files[head.source]):
            heappush(heap, vcf_files[head.source][record_pointer[head.source]])
    return ans

def ll_solve_chrom(para):
    start_time = time.time()
    vcf_files = para[0]  # vcf_files->dict{fileidx->list(Record)}
    chrom = para[1]
    anno = para[2]
    support = para[3]
    svtype = para[4]
    filenum = para[5]
    supp_filter = para[6]
    work_dir = para[7]
    config_para = para[8]
    max_dist = parse_config_para(config_para, svtype)
    sort_ans = first_sort(vcf_files, max_dist[0], filenum, (chrom, svtype))
    #print('finish sorting in %s'%(str(time.time() - start_time)))
    result = list()
    idx = 0
    for node in sort_ans:
        idx += 1
        candidates = cal_can(node, 0, svtype, max_dist)
        for candidate in candidates: # candidate -> list(Record)
            if len(candidate) < support:
                continue
            candidate_idx, cipos, ciend = cal_center(candidate)
            candidate_record = candidate[candidate_idx]  # Record
            annotation = solve_annotation(candidate_record.type, anno, candidate_record.start, candidate_record.end)  # dict{'gene_id' -> str}
            candidate_dict = dict()
            for can in candidate:
                candidate_dict[can.source] = can
            result.append([candidate_record.chrom1, candidate_record.start, candidate_record, cipos, ciend, candidate_dict, annotation])
    #print('finish %s-%s in %s, total %dnodes'%(chrom, svtype, str(time.time() - start_time), len(result)))
    output_result(result, filenum, '%stemporary/%s_%s'%(work_dir, chrom, svtype), supp_filter)
    print('finish %s-%s in %ss, total %d nodes'%(chrom, svtype, str(time.time() - start_time), len(result)))


# all variants on the chrom
def resolve_chrom(para):
    #print('start resolving chrom ' + chrom)
    vcf_files = para[0]
    chrom = para[1]
    iothread = para[4]
    start_time = time.time()
    file_dict = dict()
    type_dict = dict()
    idx = 0
    '''
    pool = Pool(processes = iothread)
    for filename in vcf_files:
        file_dict[idx] = pool.map_async(parse_vcf_chrom, [(filename, chrom, idx)])
        idx += 1
    pool.close()
    pool.join()
    '''
    for filename in vcf_files:
        file_dict[idx] = parse_vcf_chrom([filename, chrom, idx])
        idx += 1
    reading_time = time.time() - start_time
    #'''
    for fileidx in file_dict:
        #records = file_dict[fileidx].get()[0] # {svtype -> List(Record)}
        records = file_dict[fileidx] # {svtype -> List(Record)}
        for svtype in records:
            if svtype not in type_dict:
                type_dict[svtype] = dict()
            type_dict[svtype][fileidx] = records[svtype]
    parsing_time = time.time() - start_time
    for svtype in type_dict:
        para[0] = type_dict[svtype]
        para[4] = svtype
        ll_solve_chrom(para)
    print('finish chrom %s in %f, %f, %fs'%(chrom, reading_time, parsing_time, time.time() - start_time))


def conquer_merge(vcf_files, chrom, anno, support, svtype, filenum, supp_filter, work_dir, config_para):
    start_time = time.time()
    max_dist = parse_config_para(config_para, svtype)
    sort_ans = first_sort(vcf_files, max_dist[0], filenum, (chrom, svtype))
    #print('finish sorting in %s'%(str(time.time() - start_time)))
    result = list()
    idx = 0
    for node in sort_ans:
        idx += 1
        candidates = cal_can(node, 0, svtype, max_dist)
        for candidate in candidates: # candidate -> list(Record)
            if len(candidate) < support:
                continue
            candidate_idx, cipos, ciend = cal_center(candidate)
            candidate_record = candidate[candidate_idx]  # Record
            annotation = solve_annotation(candidate_record.type, anno, candidate_record.start, candidate_record.end)  # dict{'gene_id' -> str}
            candidate_dict = dict()
            for can in candidate:
                candidate_dict[can.source] = can
            result.append([candidate_record.chrom1, candidate_record.start, candidate_record, cipos, ciend, candidate_dict, annotation])
    #print('finish %s-%s in %s, total %dnodes'%(chrom, svtype, str(time.time() - start_time), len(result)))
    output_result(result, filenum, '%stemporary/%s_%s'%(work_dir, chrom, svtype), supp_filter)
    print('finish %s-%s in %ss, total %d nodes'%(chrom, svtype, str(time.time() - start_time), len(result)))

def fighting(args):
    start_time = time.time()
    if args.work_dir[-1] != '/':
        args.work_dir += '/'
    if not os.path.exists('%stemporary'%(args.work_dir)):
        os.mkdir('%stemporary'%(args.work_dir))
    else:
        print('temporary directory existed.')
        exit(0)
    if not os.path.exists('%sindex'%(args.work_dir)):
        os.mkdir('%sindex'%(args.work_dir))
    else:
        print('index directory existed.')
        exit(0)
    filenames = index_vcf(args.input, args.threads, args.work_dir)
    annotation_dict = parse_annotation_file(args.annotation)
    filter_chrom_list = None if args.filter_chrom is None else args.filter_chrom.split(',')
    if args.wgs:
        config_para = parse_para(os.path.abspath(sys.argv[0])[:(os.path.abspath(sys.argv[0]).rfind('/') + 1)] + 'config_WGS')
    else:
        config_para = parse_para(os.path.abspath(sys.argv[0])[:(os.path.abspath(sys.argv[0]).rfind('/') + 1)] + 'config')
    if args.max_dist != -1:
        for i in range(5):
            config_para[i] = args.max_dist
    
    #sample_ids, contiginfo, contigs = pre_vcfs(filenames, args.IOthreads, filter_chrom_list)
    sample_ids, contiginfo = resolve_contigs(filenames, args.IOthreads, filter_chrom_list)
    print('finish parsing in %ss'%(round(time.time() - start_time, 6)))
    file = open(args.output, 'w')
    generate_header(file, contiginfo, sample_ids)
    file.close()

    pool = Pool(processes = args.threads)
    para = []
    for contig in contiginfo:
        if annotation_dict != None and contig in annotation_dict:
            anno = annotation_dict[contig]
        else:
            anno = []
        
        file_dict = dict()
        type_dict = dict()
        idx = 0
        ss_time = time.time()
        for filename in filenames:
            file_dict[idx] = parse_vcf_chrom([filename, contig, idx])  # 这里可不可以开多进程？
            idx += 1
        #'''
        for fileidx in file_dict:
            #records = file_dict[fileidx].get()[0] # {svtype -> List(Record)}
            records = file_dict[fileidx] # {svtype -> List(Record)}
            for svtype in records:
                if svtype not in type_dict:
                    type_dict[svtype] = dict()
                type_dict[svtype][fileidx] = records[svtype]
        print(time.time() - ss_time)
        for svtype in type_dict:
            #para.append([type_dict[svtype], contig, anno, args.support, svtype, len(sample_ids), args.supp, args.work_dir, config_para])
            multi_para = [(type_dict[svtype], contig, anno, args.support, svtype, len(sample_ids), args.supp, args.work_dir, config_para)]
            pool.map_async(multi_run_wrapper, multi_para)

    #pool.map(ll_solve_chrom, para)
    pool.close()
    pool.join()
    print('finish merging in %ss'%(round(time.time() - start_time, 6)))
    os.system('cat %stemporary/* > %stemporary/vcf'%(args.work_dir, args.work_dir))
    os.system('sort -k 1,1 -k 2,2n %stemporary/vcf >> %s'%(args.work_dir, args.output))
    os.system('rm -r %sindex'%(args.work_dir))
    os.system('rm -r %stemporary'%(args.work_dir))

    #os.system('rm %stemporary/vcf'%(args.work_dir))
    print('finish in %ss'%(round(time.time() - start_time, 6)))
def multi_run_wrapper(args):
	return conquer_merge(*args)

# move io into process
def main_ctrl(args):
    start_time = time.time()
    if args.work_dir[-1] != '/':
        args.work_dir += '/'
    if not os.path.exists('%stemporary'%(args.work_dir)):
        os.mkdir('%stemporary'%(args.work_dir))
    else:
        print('temporary directory existed.')
        exit(0)
    if not os.path.exists('%sindex'%(args.work_dir)):
        os.mkdir('%sindex'%(args.work_dir))
    else:
        print('index directory existed.')
        exit(0)
    filenames = index_vcf(args.input, args.threads, args.work_dir)
    annotation_dict = parse_annotation_file(args.annotation)
    filter_chrom_list = None if args.filter_chrom is None else args.filter_chrom.split(',')
    if args.wgs:
        config_para = parse_para(os.path.abspath(sys.argv[0])[:(os.path.abspath(sys.argv[0]).rfind('/') + 1)] + 'config_WGS')
    else:
        config_para = parse_para(os.path.abspath(sys.argv[0])[:(os.path.abspath(sys.argv[0]).rfind('/') + 1)] + 'config')
    if args.max_dist != -1:
        for i in range(5):
            config_para[i] = args.max_dist
    
    sample_ids, contiginfo = resolve_contigs(filenames, args.IOthreads, filter_chrom_list)
    print('finish parsing in %ss'%(round(time.time() - start_time, 6)))
    #print('%d samples find'%(len(sample_ids)))
    file = open(args.output, 'w')
    generate_header(file, contiginfo, sample_ids)
    file.close()
    #os.system('mkdir %stemporary'%(args.work_dir))
    pool = Pool(processes = args.threads)
    para = []
    for contig in contiginfo:
        if annotation_dict != None and contig in annotation_dict:
            anno = annotation_dict[contig]
        else:
            anno = []
        #resolve_chrom([filenames, contig, anno, args.support, args.IOthreads, len(sample_ids), args.supp, args.work_dir, config_para])
        para.append([filenames, contig, anno, args.support, args.IOthreads, len(sample_ids), args.supp, args.work_dir, config_para])
    pool.map(resolve_chrom, para)
    pool.close()
    pool.join()
    print('finish merging in %ss'%(round(time.time() - start_time, 6)))
    os.system('cat %stemporary/* > %stemporary/vcf'%(args.work_dir, args.work_dir))
    os.system('sort -k 1,1 -k 2,2n %stemporary/vcf >> %s'%(args.work_dir, args.output))
    #os.system('rm -r %sindex'%(args.work_dir))
    os.system('rm -r %stemporary'%(args.work_dir))

    #os.system('rm %stemporary/vcf'%(args.work_dir))
    print('finish in %ss'%(round(time.time() - start_time, 6)))

# not seperate chrom, all read in memory
def main(args):
    start_time = time.time()
    chr_dict, sample_ids, contiginfo, order = parse_vcfs(args.input, args.IOthreads, args.filter_chrom)
    if args.supp != None and len(sample_ids) != len(args.supp):
        print('parameter \'supp\' doesn\'t match number of files input in %s'%(args.input))
        exit(0)
    if args.work_dir[-1] != '/':
        args.work_dir += '/'
    annotation_dict = parse_annotation_file(args.annotation)
    if args.wgs:
        config_para = parse_para(os.path.abspath(sys.argv[0])[:(os.path.abspath(sys.argv[0]).rfind('/') + 1)] + 'config_WGS')
    else:
        config_para = parse_para(os.path.abspath(sys.argv[0])[:(os.path.abspath(sys.argv[0]).rfind('/') + 1)] + 'config')
    if args.max_dist != -1:
        for i in range(5):
            config_para[i] = args.max_dist

    print('finish parsing in %s'%(str(time.time() - start_time)))
    file = open(args.output, 'w')
    generate_header(file, contiginfo, sample_ids)
    file.close()

    pool = Pool(processes = args.threads)
    os.system('mkdir %stemporary'%(args.work_dir))
    para = []
    for chr_o in order:
        chr = chr_o[0]
        if annotation_dict != None and chr in annotation_dict:
            anno = annotation_dict[chr]
        else:
            anno = []
        for svtype in chr_dict[chr]:
            para.append([chr_dict[chr][svtype], chr, anno, args.support, svtype, len(sample_ids), args.supp, args.work_dir, config_para])
    #para.append([chr_dict['2']['DEL'], '2', None, args.support, 'DEL', len(sample_ids), args.supp])
    pool.map(ll_solve_chrom, para)
    pool.close()
    pool.join()
    os.system('cat %stemporary/* > %stemporary/vcf'%(args.work_dir, args.work_dir))
    os.system('sort -k 1,1 -k 2,2n %stemporary/vcf >> %s'%(args.work_dir, args.output))
    #os.system('rm -r ' + args.work_dir + 'index/')
    os.system('rm -r %stemporary'%(args.work_dir))
    print('finish in ' + str(time.time() - start_time) + 'seconds')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
            metavar = 'input_file',
            type = str,
            help = 'filenames of VCF files to be merged')
    parser.add_argument('output',
            type = str,
            help = 'output VCF format file')
    parser.add_argument('work_dir',
            type = str,
            help = 'work directory for temporary files')
    parser.add_argument('-t', '--threads',
            type = int,
            default = 16,
            help = 'number of threads to use[%(default)s]')
    parser.add_argument('-tio', '--IOthreads',
            type = int,
            default = 2,
            help = 'number of threads to use[%(default)s]')
    parser.add_argument('-a', '--annotation',
            type = str,
            default = None,
            help = 'annotation file to add')
    parser.add_argument('--massive',
            action="store_true",
            default = False)
    parser.add_argument('-S', '--support',
            type = int,
            default = 1,
            help = 'support vector number[%(default)s]')
    parser.add_argument('-max_dist',
            type = int,
            default = -1,
            help = 'Maximum distance[%(default)s]')
    parser.add_argument('-max_ins_ratio',
            type = float,
            default = 0.3,
            help = 'Maximum ratio of merging insertions[%(default)s]')
    parser.add_argument('--comp_strand',
            action="store_true",
            default = False)
    parser.add_argument('-supp',
            type = str,
            default = None,
            help = 'support vector to filter')
    parser.add_argument('-filter_chrom',
            type = str,
            default = None,
            help = 'only handle these chromosomes, separate by comma and no space')
    parser.add_argument('--wgs',
            action="store_true",
            default = False,
            help = 'merging SVs from WGS callers')
 
    args = parser.parse_args(sys.argv[1:])
    if args.massive:
        main_ctrl(args)
        #fighting(args)
    else:
        main(args)
    # /usr/bin/time -v python first_release/cuteSV_merge.py real_test/real_100.fofn temp ./

'''
chrom_set = set()
base_cmd = 'cat '
for vcf_filename in vcf_filenames:
    base_cmd += vcf_filename + ' '
base_cmd += '| grep -v \'#\' | awk -F \'\\t\' \'{print $1}\' | uniq | sort | uniq'
print(base_cmd)
'''
'''
chrom_set = set()
contiginfo = dict()
sample_set = set()
sample_ids = list()
for vcf_filename in vcfgz_filenames:  # 默认header的config中存放了所有chrom信息
    vcf_reader = VariantFile(vcf_filename, 'r')
    for i in range(len(vcf_reader.header.contigs)):
        chrom_set.add(str(vcf_reader.header.contigs[i].name))
        contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
    sample_id = vcf_reader.header.samples[0]
    temp_sample_id = sample_id
    temp_idx = 0
    while temp_sample_id in sample_set:
        temp_sample_id = sample_id + '_' + str(temp_idx)
        temp_idx += 1
    sample_ids.append(temp_sample_id)
    sample_set.add(temp_sample_id)
chrom_cnt = []

base_cmd1 = 'grep -v \'#\' ' + vcf_filenames[0] + ' | awk -F \'\\t\' \'{print $1}\' | grep -x \''
base_cmd2 = '\' | wc -l'
for chrom in chrom_set:
    fd = os.popen(base_cmd1 + chrom + base_cmd2)
    mi = fd.read().strip()
    chrom_cnt.append([chrom, int(mi)])
chrom_cnt.sort(key = lambda x:x[1], reverse = True)
for chrom in chrom_set:
    chrom_cnt.append([chrom])
'''
'''
for iter in chrom_cnt:
    if annotation_dict != None and iter[0] in annotation_dict:
        anno = annotation_dict[iter[0]]
    else:
        anno = []
    result.append(ll_solve_chrom([vcfgz_filenames, iter[0], max_dist, max_inspro, anno]))
'''
'''
for iter in chrom_cnt:
    # multi processes
    #print(iter[0])
    #if iter[1] == 0:
    #    continue
    if annotation_dict != None and iter[0] in annotation_dict:
        anno = annotation_dict[iter[0]]
    else:
        anno = []
    para = [(vcfgz_filenames, iter[0], max_dist, max_inspro, anno)]
    result.append(pool.map_async(ll_solve_chrom, para))
'''
'''
para = []
for iter in chrom_cnt:
    if annotation_dict != None and iter[0] in annotation_dict:
        anno = annotation_dict[iter[0]]
    else:
        anno = []
    para.append([vcfgz_filenames, iter[0], max_dist, max_inspro, anno])
result = pool.map(ll_solve_chrom, para, chunksize = 1)
'''
'''
def avl_solve_chrom(vcf_filenames, chrom):
    tree = AVLTree()
    #print('start solving chrom: ' + chrom)
    time0 = time.time()
    for i in range(len(vcf_filenames)):
        idx = 0
        vcf_reader = VariantFile(vcf_filenames[i], 'r')
        for record in vcf_reader.fetch(chrom):
            idx += 1
            tree.insert(i, Record(record, i))
    time1 = time.time()
    output_chrom(tree, len(vcf_filenames), chrom)
    time2 = time.time()
'''
