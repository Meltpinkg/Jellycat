run_sim() {
    #run_sim dir
    for i in {5x,10x,15x,20x,30x}
	do
		python first_release/cuteSV_merge.py input/cutesv_$i.fofn $1/sim_$i.vcf ./ -t 16
        bgzip -c $1/sim_$i.vcf > $1/sim_$i.vcf.gz
        tabix $1/sim_$i.vcf.gz
        rm -r cmp_$i
        nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/sim_$i.vcf.gz -o cmp_$i -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	done
}

run_test() {
    rm -r $1
    mkdir $1
    rm answer.txt
    echo 'simulated:'
    bash first_release/benchmark/run_merge_sample.sh $1
    #run_sim $1
    #python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt answer.txt
    echo 'TRIO:'
    python first_release/cuteSV_merge.py trio/trio_ccs.fofn $1/trio_mer.vcf ./ -t 16
    python first_release/benchmark/test_trio.py $1/trio_mer.vcf $1/trio_mer.tmp answer.txt
    echo 'HG002:'
    python first_release/cuteSV_merge.py HG002/HG002_cuteSV.fofn $1/2_mer.vcf ./ -t 16
    python first_release/benchmark/test_HG002.py $1/2_mer.vcf answer.txt
    rm nohup.out
}
:<<!
for it in {0.3,0.4,0.6,0.7}
do
    for j in {0.3,0.5,0.6,0.8}
    do
        for k in {0.3,0.4,0.5,0.7}
        do
            for l in {0.3,0.5,0.6,0.8}
            do
                echo $it $j $k $l >> $1
                for ii in {0,1,2,3,4}
                do
                    run_test 104 $it $j $k $l $ii $1
                done
            done
        done
    done
done
!
run_test 0129

run_jas() {
    for i in {5x,10x,15x,20x,30x}
    do
        jasmine file_list=input/cutesv_$i.fofn out_file=output/jasmine_cutesv_$i.vcf max_dist_linear=0.5 min_dist=1000
        bash src/sort_vcf.sh output/jasmine_cutesv_$i.vcf
        bgzip -c output/jasmine_cutesv_$i.vcf > output/jasmine_cutesv_$i.vcf.gz
        tabix output/jasmine_cutesv_$i.vcf.gz
        rm -r cmp_$i
        nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c output/jasmine_cutesv_$i.vcf.gz -o cmp_$i -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
    done
    rm nohup.out
    python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt answer.txt
}
run_sur() {
    for i in {5x,10x,15x,20x,30x}
    do
        /data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/cutesv_$i.fofn 1000 2 1 0 1 30 output/survivor_cutesv_$i.vcf
        bash src/sort_vcf.sh output/survivor_cutesv_$i.vcf
        bgzip -c output/survivor_cutesv_$i.vcf > output/survivor_cutesv_$i.vcf.gz
        tabix output/survivor_cutesv_$i.vcf.gz
        rm -r cmp_$i
        nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c output/survivor_cutesv_$i.vcf.gz -o cmp_$i -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
    done
    rm nohup.out
    python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt answer.txt
}
