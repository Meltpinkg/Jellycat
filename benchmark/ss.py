import sys
import json
#print(sys.argv[1])
#print(sys.argv[2])


length = len(sys.argv)
file = open(sys.argv[length - 1], 'a')
mat = []
for i in range(1, length - 1) :
    with open(sys.argv[i], 'r') as f:
        dct = json.load(f)
        anslist = [dct['precision'], dct['recall'], dct['f1'], dct['TP-base'], dct['TP-call'], dct['TP-base']+dct['FN'], dct['TP-call']+dct['FP']] # precision recall F1 TP-base TP-call TP-base+FN TP-call+FP
        
        for x in range(0, 7):
            file.write(str(anslist[x])+',')
        file.write('\n')
        
        mat.append([])
        for x in range(0, 3):
            #testtxt.write(str(anslist[x])+',')
            mat[-1].append(anslist[x])
file.write('\n')
#print(mat)
#print('%f,%f,%f,%f,%f,%f,%f'%(mat[0][2],mat[1][2],mat[2][2],mat[3][2],mat[4][0],mat[4][1],mat[4][2]))
        #testtxt.write('\n')
