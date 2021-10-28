import sys

#print(sys.argv[1])
#print(sys.argv[2])


length = len(sys.argv)
file = open(sys.argv[length - 1], 'a')
mat = []
for i in range(1, length - 1) :
    with open(sys.argv[i], 'r') as f:
        anslist = list()
        anslist = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # anslist_gt = list()
        # anslist_gt = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        for line in f:
            call_number = 0
            base_number = 0
            if 'precision' in line and 'gt_precision' not in line:
                a = line.split(':')[1].split(',')[0]
                try:
                    anslist[0] = float(a)
                except:
                    anslist[0] = 0
            if 'recall' in line and 'gt_recall' not in line:
                a = line.split(':')[1].split(',')[0]
                try:
                    anslist[1] = float(a)
                except:
                    anslist[1] = 0
            if 'f1' in line and 'gt_f1' not in line:
                a = line.split(':')[1].split(',')[0]
                try:
                    anslist[2] = float(a)
                except:
                    anslist[2] = 0
            if 'TP-base' in line and 'TP-base_TP-gt' not in line and "TP-base_FP-gt" not in line:
                a = line.split(':')[1].split(',')[0]
                anslist[3] = int(a)
            if 'TP-call' in line and 'TP-call_TP-gt' not in line and "TP-call_FP-gt" not in line:
                a = line.split(':')[1].split(',')[0]
                anslist[4] = int(a)
            if 'call cnt' in line:
                a = line.split(':')[1].split(',')[0]
                anslist[5] = int(a)
            if 'base cnt' in line:
                a = line.split(':')[1].split(',')[0]
                anslist[6] = int(a)

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
