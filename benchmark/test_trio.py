from pysam import VariantFile
import sys
lst = [90272, 139408, 533483, 1219553, 1224402, 2054692, 2256563, 2583572, 2594415, 2596458, 2604982, 2613655, 2623286, 2634216, 2684816, 2954350, 3344757, 3393869, 3393869, 3561912, 7580282, 9684000, 17124911, 18376200, 19386826, 19387909, 19387914, 23603557, 26534275, 28399624, 29691912, 31904977, 35267851, 36124478, 36679545, 38405305, 39677664, 40130171, 45688745, 46537404]

#gt: record.samples[i]['GT'] (tuple)
def parse_gt(gt):
    if gt == None or gt[0] == None:
        return 0
    if gt[0] + gt[1] == 0:
        return 0
    if gt[0] + gt[1] == 1:
        return 1
    if gt[0] + gt[1] == 2:
        return 2

# gt:  1/1:A:ACGGGTGCC
def parse_gt_inline(gt):
    cnt = 0
    if gt[0] == '1':
        cnt += 1
    if gt[2] == '1':
        cnt += 1
    return cnt

# False -> inconsistency
def check(son, fa, ma):
    if fa == 0 and ma == 0:
        if son == 1 or son == 2:
            return False
    if fa == 0 and ma == 1 and son == 2:
        return False
    if fa == 1 and ma == 0 and son == 2:
        return False
    if fa == 0 and ma == 2:
        if son == 0 or son == 2:
            return False
    if fa == 2 and ma == 0:
        if son == 0 or son == 2:
            return False
    if fa == 1 and ma == 2 and son == 0:
        return False
    if fa == 2 and ma == 1 and son == 0:
        return False
    if fa == 2 and ma == 2:
        if son == 0 or son == 1:
            return False
    return True

def main(vcf_file, output_file):
    vcf_reader = VariantFile(vcf_file, 'r')
    bcf_out = VariantFile(output_file, 'w', header=vcf_reader.header)
    '''
    sample1 = vcf_reader.header.samples[0]
    sample2 = vcf_reader.header.samples[1]
    sample3 = vcf_reader.header.samples[2]
    '''
    cnt = 0
    tot = 0
    for record in vcf_reader.fetch():
        tot += 1
        #supp_vec = record.info['SUPP_VEC']
        son = parse_gt(record.samples[0]['GT'])
        father = parse_gt(record.samples[1]['GT'])
        mother = parse_gt(record.samples[2]['GT'])
        res = check(son, father, mother)
        if res == False:
            cnt += 1
            bcf_out.write(record)
    print(cnt)
    print(tot)
    print(cnt / tot)
    return cnt, tot

def main_indel(vcf_file, output_file):
    output = open(output_file, 'w')
    testtxt = open('answer.txt', 'a')
    cnt = dict()
    tot = dict()
    with open(vcf_file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            seq = line.strip().split('\t')
            son = parse_gt_inline(seq[-3])
            father = parse_gt_inline(seq[-2])
            mother = parse_gt_inline(seq[-1])
            svtype = line.split('SVTYPE=')[1].split(';')[0]
            if svtype not in tot:
                tot[svtype] = 0
            tot[svtype] += 1
            res = check(son, father, mother)
            if res == False:
                if svtype not in cnt:
                    cnt[svtype] = 0
                cnt[svtype] += 1
                output.write('=====')
            output.write(line)
    ccnt = 0
    ttot = 0
    for svtype in ['INS', 'DEL', 'INV', 'DUP', 'BND']:
        if svtype == 'BND' and 'BND' not in cnt:
            svtype = 'TRA'
        print('%s %f'%(svtype, cnt[svtype] / tot[svtype]))
        testtxt.write('%f,'%(cnt[svtype] / tot[svtype]))
        ccnt += cnt[svtype]
        ttot += tot[svtype]
    if ttot == 0:
        rat = -1
    else:
        rat = ccnt/ttot
    testtxt.write('%d,%d,%f\n'%(ccnt, ttot, rat))
    print('%d,%d,%f'%(ccnt, ttot, rat))
    testtxt.close()

if __name__ == '__main__':
    #main(sys.argv[1], sys.argv[2])
    main_indel(sys.argv[1], sys.argv[2])