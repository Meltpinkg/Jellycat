from pysam import VariantFile
import sys

def main(vcf_file):
    testtxt = open('answer.txt', 'a')
    supp_dict = dict()
    vcf_reader = VariantFile(vcf_file, 'r')
    tot = 0
    for record in vcf_reader.fetch():
        supp_vec = record.info['SUPP_VEC']
        tot += 1
        if supp_vec in supp_dict:
            supp_dict[supp_vec] = supp_dict[supp_vec] + 1
        else:
            supp_dict[supp_vec] = 1
    #print('test HG002')
    for i in ['001', '011', '101', '111', '110', '010', '100']:
        print(i + ': ' + str(supp_dict[i]))
        testtxt.write('%d,'%(supp_dict[i]))
    '''
    print(supp_dict['001'])
    print(supp_dict['011'])
    print(supp_dict['101'])
    print(supp_dict['111'])
    print(supp_dict['110'])
    print(supp_dict['010'])
    print(supp_dict['100'])
    '''
    print(supp_dict['011'] + supp_dict['101'] + supp_dict['110'] + supp_dict['111'])
    print(tot)
    print((supp_dict['011'] + supp_dict['101'] + supp_dict['110'] + supp_dict['111']) / (tot))
    calmdown = supp_dict['011'] + supp_dict['101'] + supp_dict['110'] + supp_dict['111']
    testtxt.write('%d,%f\n'%(calmdown, calmdown/tot))
    testtxt.close()

if __name__ == '__main__':
    main(sys.argv[1])