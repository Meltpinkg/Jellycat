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
#from numpy.core.fromnumeric import sort
from cuteSV_calculation import *
from cuteSV_linkedList import Record, parse_svtype
from cuteSV_output import output_result, solve_annotation, generate_header
from pysam import VariantFile
from multiprocessing import Pool
from heapq import *
import os
import sys
import time
import argparse

def create_index(line, gz_filename):
    try:
        cmd = 'bgzip -c ' + line + ' > ' + gz_filename
        os.system(cmd)
        os.system('tabix -p vcf ' + gz_filename)
    except Exception as ee:
        print(ee)

# index vcf files and output .gz in work_dir
def index_vcf(filenames, threads, work_dir):
    #print('start indexing...')
    start_time = time.time()
    if work_dir[-1] != '/':
        work_dir += '/'
    if not os.path.exists(work_dir + 'index'):
        os.mkdir(work_dir + 'index')
    else:
        print('Index directory existed.')
    process_pool = Pool(processes = threads)
    vcf_filenames = []
    vcfgz_filenames = []
    with open(filenames, 'r') as f:
        for line in f:
            line = line.strip()
            vcf_filenames.append(line)
            filename = line.split('/')[-1]
            if filename[-3:] == 'vcf':
                filename += '.gz'
                gz_filename = work_dir + line.split('/')[-1] + '.gz'
                vcfgz_filenames.append(gz_filename)
                process_pool.apply_async(create_index, (line, gz_filename))
            elif filename[-6:] == 'vcf.gz':
                vcfgz_filenames.append(line)
            else:
                print('input error')
                continue
            
    process_pool.close()
    process_pool.join()
    print('finish indexing in %s'%(round(time.time() - start_time, 6)))
    #return vcfgz_filenames, chrom_cnt, contiginfo, sample_ids
    return vcfgz_filenames

# parse vcf and return {chr -> list(Record)}
def parse_vcf(para):
    filename = para[0]
    idx = para[1]
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict()
    # contigs
    contiginfo = dict()
    for i in range(len(vcf_reader.header.contigs)):
        try:
            contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[str(vcf_reader.header.contigs[i].name)] = 0
    record_dict['contig'] = contiginfo
    # records
    for record in vcf_reader.fetch():
        chr = record.chrom
        svtype = parse_svtype(record.info['SVTYPE'])
        if chr not in record_dict:
            record_dict[chr] = dict()
        if svtype not in record_dict[chr]:
            record_dict[chr][svtype] = []
        record_dict[chr][svtype].append(Record(record, idx))
    # sample_id
    record_dict['sample'] = vcf_reader.header.samples[0]
    return record_dict

# parse vcfs and return {chr -> {fileidx -> List(Record)}}
def parse_vcfs(filenames, threads):
    print('start parse vcfs...')
    pool = Pool(processes = threads)
    file_dict = dict()
    idx = 0
    with open(filenames, 'r') as f:
        for filename in f:
            file_dict[idx] = pool.map_async(parse_vcf, [(filename.strip(), idx)])
            idx += 1
    pool.close()
    pool.join()

    chr_dict = dict()
    chr_cnt = dict()
    sample_temp = ['' for i in range(idx)]
    contiginfo = dict()
    for fileidx in file_dict:
        records = file_dict[fileidx].get()[0]  #{chr -> svtype -> List(Record)}
        sample_temp[fileidx] = records['sample']
        contig_temp = records['contig']
        records.pop('sample')
        records.pop('contig')
        for chr in records:
            if chr not in chr_dict:
                chr_dict[chr] = dict()
                chr_cnt[chr] = 0
            if chr not in contiginfo:
                contiginfo[chr] = contig_temp[chr]
            for svtype in records[chr]:
                if svtype not in chr_dict[chr]:
                    chr_dict[chr][svtype] = dict()
                chr_dict[chr][svtype][fileidx] = records[chr][svtype]
                chr_cnt[chr] += len(records[chr][svtype])
    order = sorted(chr_cnt.items(), key=lambda x:x[1], reverse=True) # List
    # rename samples
    sample_set = set()
    sample_ids = list()
    for sample_id in sample_temp:
        temp_sample_id = sample_id
        temp_idx = 0
        while temp_sample_id in sample_set:
            temp_sample_id = sample_id + '_' + str(temp_idx)
            temp_idx += 1
        sample_ids.append(temp_sample_id)
        sample_set.add(temp_sample_id)
    return chr_dict, sample_ids, contiginfo, order

# all variants on the chrom
def resolve_chrom(filenames, output, chrom, sample_ids, threads, max_dist, max_inspro, seperate_svtype, anno, support, diff_ratio_merging_INS, diff_ratio_merging_DEL):
    #print('start resolving chrom ' + chrom)
    start_time = time.time()
    pool = Pool(processes = threads)
    record_dict = dict()
    idx = 0
    for filename in filenames:
        record_dict[idx] = pool.map_async(parse_vcf_chrom, [(filename, chrom, idx)])
        idx += 1
    pool.close()
    pool.join()
    '''
    if seperate_svtype:
        result = solve_chrom2(record_dict, max_dist, max_inspro, anno, support, diff_ratio_merging_INS, diff_ratio_merging_DEL)
    else:
        result = solve_chrom1(record_dict, max_dist, max_inspro, anno, support, diff_ratio_merging_INS, diff_ratio_merging_DEL)
    result = sorted(result, key = lambda x:(x[0], int(x[1])))
    output_result(result, len(sample_ids), output)
    print('finish resolving %s, %d in %s'%(chrom, len(result), round(time.time() - start_time, 4)))
    '''

def parse_vcf_chrom(para):
    filename = para[0]
    chrom = para[1]
    idx = para[2]
    vcf_reader = VariantFile(filename, 'r')
    record_list = list()
    # records
    for record in vcf_reader.fetch(chrom):
        record_list.append(Record(record, idx))
    return record_list

def resolve_contigs(filenames, threads):
    print('start resolving contigs...')
    start_time = time.time()
    sample_temp = list()
    sample_set = set()
    sample_ids = list()
    contiginfo = dict()
    result = list()
    pool = Pool(processes = threads)
    for filename in filenames:
        result.append(pool.map_async(parse_contigs, [(filename)]))
    pool.close()
    pool.join()
    for res in result:
        temp = res.get()[0]
        for contig in temp:
            if contig == 'sample':
                sample_id = temp['sample']
                temp_sample_id = sample_id
                temp_idx = 0
                while temp_sample_id in sample_set:
                    temp_sample_id = sample_id + '_' + str(temp_idx)
                    temp_idx += 1
                sample_ids.append(temp_sample_id)
                sample_set.add(temp_sample_id)
            else:
                if contig not in contiginfo:
                    contiginfo[contig] = temp[contig]
    print('finish resolving contigs in %s'%(round(time.time() - start_time, 4)))
    return sample_ids, contiginfo

def parse_contigs(para):
    filename = para
    contiginfo = dict()
    vcf_reader = VariantFile(filename, 'r')
    for i in range(len(vcf_reader.header.contigs)):
        contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
    contiginfo['sample'] = vcf_reader.header.samples[0]
    return contiginfo

def parse_annotation_file(annotation_file):
    if annotation_file == None:
        return None
    annotation_dict = dict()  # [chrom -> [[a, b, dict()], ...]]
    start_time = time.time()
    with open(annotation_file, 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            chr = seq[0]
            if chr[:3] == "chr":
                chr = chr[3:]
            if chr not in annotation_dict:
                annotation_dict[chr] = []
            info_seq = seq[8].strip().split(';')
            info_dict = dict()
            for info_item in info_seq[:-1]:
                info_item = info_item.strip().split(' ')
                if len(info_item) != 2:
                    print('error:length of gtf_info_item is not 2')
                    continue
                info_dict[info_item[0]] = info_item[1][1:-1]
            annotation_dict[chr].append([int(seq[3]), int(seq[4]), info_dict])    
    print(time.time() - start_time)
    return annotation_dict

# vcf_files->dict{fileidx->list(Record)}  only one svtype
# return list[]
def first_sort(vcf_files, max_dist, filenum):
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
    vcf_files = para[0]
    chrom = para[1]
    anno = para[2]
    support = para[3]
    svtype = para[4]
    filenum = para[5]
    max_dist = 1000
    sort_ans = first_sort(vcf_files, max_dist, filenum)
    #print('finish sorting in %s'%(str(time.time() - start_time)))
    result = list()
    idx = 0
    for node in sort_ans:
        idx += 1
        if svtype == 'INS':
            candidates = cal_can_INS(node, 0)
        elif svtype == 'DEL':
            candidates = cal_can_DEL(node, 0)
        elif svtype == 'INV':
            candidates = cal_can_INV(node, 0)
        elif svtype == 'DUP':
            candidates = cal_can_DUP(node, 0)
        else:
            candidates = cal_can_BND(node, 0)
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
    print('finish %s-%s in %s'%(chrom, svtype, str(time.time() - start_time)))
    output_result(result, filenum, 'temporary/' + chrom + '_' + svtype)
    print('finish %s-%s in %s, total %dnodes'%(chrom, svtype, str(time.time() - start_time), len(result)))

# resolve_chrom
def main_ctrl(args):
    start_time = time.time()
    max_dist = args.max_dist
    max_inspro = args.max_inspro
    filenames = index_vcf(args.input, args.threads, args.work_dir + 'index/')
    annotation_dict = parse_annotation_file(args.annotation)
    sample_ids, contiginfo = resolve_contigs(filenames, args.IOthreads)
    print('%d samples find'%(len(sample_ids)))
    file = open(args.output, 'w')
    generate_header(file, contiginfo, sample_ids)
    file.close()
    for contig in contiginfo:
        if annotation_dict != None and contig in annotation_dict:
            anno = annotation_dict[contig]
        else:
            anno = []
        resolve_chrom(filenames, args.output, contig, sample_ids, args.IOthreads, max_dist, max_inspro, args.seperate_svtype, anno, args.support, args.diff_ratio_merging_INS, args.diff_ratio_merging_DEL)
    print('finish merging in %s'%(round(time.time() - start_time, 4)))
    os.system('rm -r ' + args.work_dir + 'index/')

# not seperate chrom, all read in memory
def main(args):
    start_time = time.time()
    chr_dict, sample_ids, contiginfo, order = parse_vcfs(args.input, args.IOthreads)
    annotation_dict = parse_annotation_file(args.annotation)
    #max_dist, max_ratio = parse_para(args)
    print('finish parsing in %s'%(str(time.time() - start_time)))
    file = open(args.output, 'w')
    generate_header(file, contiginfo, sample_ids)
    file.close()

    pool = Pool(processes = args.threads)
    os.system('mkdir temporary')
    para = []
    for chr_o in order:
        chr = chr_o[0]
        if annotation_dict != None and chr in annotation_dict:
            anno = annotation_dict[chr]
        else:
            anno = []
        for svtype in chr_dict[chr]:
            para.append([chr_dict[chr][svtype], chr, anno, args.support, svtype, len(sample_ids)])
    pool.map(ll_solve_chrom, para)
    pool.close()
    pool.join()
    os.system('cat temporary/* > temporary/vcf')
    os.system('sort -k 1,1 -k 2,2n temporary/vcf >> ' + args.output)
    #os.system('rm -r ' + args.work_dir + 'index/')
    os.system('rm -r temporary')
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
            default = 16,
            help = 'number of threads of IO[%(default)s]')
    parser.add_argument('-a', '--annotation',
            type = str,
            default = None,
            help = 'annotation file to add')
    parser.add_argument('--seperate_svtype',
            action="store_true",
            default = False)
    parser.add_argument('--seperate_chrom',
            action="store_true",
            default = False)
    parser.add_argument('-S', '--support',
            type = int,
            default = 1,
            help = 'support vector number[%(default)s]')
    parser.add_argument('-max_ins_dist',
            type = int,
            default = 1000,
            help = 'Maximum distance[%(default)s]')
    parser.add_argument('-max_ins_ratio',
            type = float,
            default = 0.3,
            help = 'Maximum ratio of merging insertions[%(default)s]')
 
    args = parser.parse_args(sys.argv[1:])
    if args.seperate_chrom:
        main_ctrl(args)
    else:
        main(args)
    # /usr/bin/time -v python src/cuteSV_merge.py input/file_SRR.fofn sample_merged25.vcf ./ -t 16 -a hg19.refGene.gtf
    # /usr/bin/time -v python src/cuteSV_merge.py test.info test_merged.vcf ./ -t 16 -a hg19.refGene.gtf
    # /usr/bin/time -v python src/cuteSV_merge.py input/cutesv.fofn merged_cutesv.vcf ./ -t 16
    # /usr/bin/time -v python src/cuteSV_merge.py input/sniffles.fofn merged1_sniffles.vcf ./ -t 16
    # /usr/bin/time -v python src/cuteSV_merge.py input/pbsv.fofn output/merged_pbsv.vcf ./ -t 16

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
