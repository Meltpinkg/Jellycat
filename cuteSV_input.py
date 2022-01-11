from cuteSV_linkedList import Record, parse_svtype
from multiprocessing import Pool
from pysam import VariantFile
import os
import time

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
    process_pool = Pool(processes = threads)
    vcfgz_filenames = []
    with open(filenames, 'r') as f:
        for line in f:
            line = line.strip()
            filename = line.split('/')[-1]
            if filename[-3:] == 'vcf':
                filename += '.gz'
                gz_filename = work_dir + 'index/' + filename
                vcfgz_filenames.append(gz_filename)
                process_pool.apply_async(create_index, (line, gz_filename)) # line -> gz_filename
            elif filename[-6:] == 'vcf.gz':
                vcfgz_filenames.append(line)
            else:
                print('input error')
                continue
    process_pool.close()
    process_pool.join()
    print('finish indexing in %ss'%(round(time.time() - start_time, 6)))
    #print(vcfgz_filenames)
    return vcfgz_filenames

# parse vcf and return {chr -> list(Record)}
def parse_vcf(para):
    filename = para[0]
    idx = para[1]
    filter_chrom_list = para[2]
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict()
    # contigs
    contiginfo = dict()
    for i in range(len(vcf_reader.header.contigs)):
        chrom = str(vcf_reader.header.contigs[i].name)
        if filter_chrom_list is not None and chrom not in filter_chrom_list:
            continue
        try:
            contiginfo[chrom] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[chrom] = 0
    record_dict['contig'] = contiginfo
    # records
    for record in vcf_reader.fetch():
        chr = record.chrom
        if filter_chrom_list is not None and chr not in filter_chrom_list:
            continue
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
def parse_vcfs(filenames, threads, filter_chrom):
    print('start parse vcfs...')
    filter_chrom_list = None
    if filter_chrom is not None:
        filter_chrom_list = filter_chrom.split(',')
    pool = Pool(processes = threads)
    file_dict = dict()
    idx = 0
    with open(filenames, 'r') as f:
        for filename in f:
            file_dict[idx] = pool.map_async(parse_vcf, [(filename.strip(), idx, filter_chrom_list)])
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

def parse_vcf_chrom(para):
    filename = para[0]
    chrom = para[1]
    idx = para[2]
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict() # svtype -> List(Record)
    # records
    try:
        for record in vcf_reader.fetch(chrom):
            try:
                svtype = parse_svtype(record.info['SVTYPE'])
            except:
                print('Warning: passing invalid SV at ' + str(record.pos))
                continue
            if svtype not in record_dict:
                record_dict[svtype] = []
            record_dict[svtype].append(Record(record, idx))
    except:
        pass
    return record_dict

def resolve_contigs(filenames, threads, filter_chrom_list):
    #print('start resolving contig infos...')
    start_time = time.time()
    sample_set = set()
    sample_ids = list()
    contiginfo = dict()
    result = list()
    pool = Pool(processes = threads)
    for filename in filenames:
        result.append(pool.map_async(parse_contigs, [(filename, filter_chrom_list)])) #file_dict[idx] = pool.map_async(parse_vcf, [(filename.strip(), idx, filter_chrom_list)])
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
    print('finish resolving contigs in %ss'%(round(time.time() - start_time, 6)))
    return sample_ids, contiginfo

def parse_contigs(para):
    filename = para[0]
    filter_chrom_list = para[1]
    contiginfo = dict()
    vcf_reader = VariantFile(filename, 'r')
    for i in range(len(vcf_reader.header.contigs)):
        chrom = str(vcf_reader.header.contigs[i].name)
        if filter_chrom_list is not None and chrom not in filter_chrom_list:
            continue
        try:
            contiginfo[chrom] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[chrom] = 0
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

def pre_vcfs(filenames, threads, filter_chrom_list):
    result = list()
    pool = Pool(processes = threads)
    for filename in filenames:
        result.append(pool.map_async(pre_vcf, [(filename, filter_chrom_list)]))
    pool.close()
    pool.join()
    re = list()
    for res in result:
        temp = res.get()[0]
        re.append(temp)
    print(re)

def pre_vcf(para):
    filename = para[0]
    filter_chrom_list = para[1]
    contiginfo = dict()
    contigs = set()
    vcf_reader = VariantFile(filename, 'r')
    for i in range(len(vcf_reader.header.contigs)):
        chrom = str(vcf_reader.header.contigs[i].name)
        if filter_chrom_list is not None and chrom not in filter_chrom_list:
            continue
        try:
            contiginfo[chrom] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[chrom] = 0
    contiginfo['sample'] = vcf_reader.header.samples[0]
    for record in vcf_reader.fetch():
        contigs.add(record.chrom)
    return contiginfo, contigs

#pre_vcfs(['/data/2/sqcao/data/merge_data/20x/CHM1_cuteSV.vcf', '/data/2/sqcao/data/merge_data/20x/CHM13_cuteSV.vcf', '/data/2/sqcao/data/merge_data/trio/HG002_CCS_svim.vcf'], 3, None)