#1 extract pccr sequences from "dm6.fa" genome by coordinates given in "panhandles.bed" file
#2 build structures of duplexes by rnaduplex
#3 extract corresponding subalignments from genome alignments (chr*.maf) of several insects
#4 build stockholm alignments of pccr regions from #3 with rnaduplex structure from #2 as annotation
from sys import argv
import subprocess
from Bio import SeqIO
from Bio.AlignIO import MafIO
import pandas as pd
import numpy as np
from ete3 import Tree
import time


#index maf files for easy search and write to dict
idx2L=MafIO.MafIndex("./maf/chr2L.mafindex", "./maf/chr2L.maf", "dm6.chr2L")
idx2R=MafIO.MafIndex("./maf/chr2R.mafindex", "./maf/chr2R.maf", "dm6.chr2R")
idx3L=MafIO.MafIndex("./maf/chr3L.mafindex", "./maf/chr3L.maf", "dm6.chr3L")
idx3R=MafIO.MafIndex("./maf/chr3R.mafindex", "./maf/chr3R.maf", "dm6.chr3R")
idx4 = MafIO.MafIndex("./maf/chr4.mafindex", "./maf/chr4.maf", "dm6.chr4")
idxM = MafIO.MafIndex("./maf/chrM.mafindex", "./maf/chrM.maf", "dm6.chrM")
idxX = MafIO.MafIndex("./maf/chrX.mafindex", "./maf/chrX.maf", "dm6.chrX")
idxY = MafIO.MafIndex("./maf/chrY.mafindex", "./maf/chrY.maf", "dm6.chrY")
ixdict = {'chr2L':idx2L,
          'chr2R':idx2R,
          'chr3L':idx3L,
          'chr3R':idx3R,
          'chr4':idx4,
          'chrM':idxM,
          'chrX':idxX,
          'chrY':idxY}


start_time_total = time.time()
genome = "./dm6.fa"   #wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
bed = argv[1]         #/home/tools/ucsc/linux.x86_64.v369/bigBedToBed panhandles_preprocessed_filtered_energy_cutoff_-15.bb panhandles.bed
tree = Tree("./dm6.27way.nh")
out = argv[2]        #where stockholm alignments will be written
out_rscape = argv[3] #output file for R-scape (unfortunately, all files i.e *.cov, *.sorted.cov and *.power together)
delete_file = argv[4] #path to outfile for R-scape (we don't use it, but still need to specify otherwize R-scape produces dozens of files)
#outdir = argv[3]  #directory for R-scape output

with open(bed, "r") as bed, open(genome, "r") as genome, open(out, "w") as out, open(out_rscape, "w") as out_rscape:
    chrs = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    iteration = 0
    start_time = time.time()
    while True:
        line = bed.readline().rstrip().split()
        if not line:
            break
        iteration+=1
#read coordinates from bed and cut pccr sequences 
        chrn = line[0]
        chromosome = chrs[chrn].seq
        pccrid = line[3]
        strand = line[5]
        length = line[10].split(",")
        
        start1 = int(line[1])
        end1 = start1 + int(length[0])
        first = chromosome[start1: end1]
        
        start2 = start1 + int(line[11].split(",")[1])
        end2 = int(line[2])
        second = chromosome[start2: end2]
        if strand == "-":
            #reverse complement
            first = first.reverse_complement()
            second = second.reverse_complement()

        header = pccrid+","+chrn+":"+str(start1)+"-"+str(end1)+"_"+str(start2)+"-"+str(end2)

#build duplex
        duplexinp = str(first) + "\\n" + str(second)
        process = subprocess.Popen("echo -e '"+duplexinp+"' | /home/tools/ViennaRNA/bin/RNAduplex",
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = process.communicate()
        duplex_header = str(output.decode())
        duplex = duplex_header.split()
#add missing points in duplex to fill in whole reference pccr sequence length
        #dc - slices of pccr handles in duplex
        dc1 = [int(i) for i in duplex[1].split(",")]
        dc2 = [int(i) for i in duplex[3].split(",")]
        duplex = duplex[0].split('&')
        duplex_left = "."*(dc1[0]-1) + duplex[0] + "."*(len(first)-dc1[1])
        duplex_right = "."*(dc2[0]-1) + duplex[1] + "."*(len(second)-dc2[1])
        
#cut subalignments from chr*.maf files
        st = 1 if strand == "+" else -1
        results1 = ixdict[chrn].get_spliced([start1], [end1], strand = st)
        results2 = ixdict[chrn].get_spliced([start2], [end2], strand = st)
        list1 = []
        list2 = []
#extract reference sequence with gaps to add corresponding gaps in duplex structure
        for seq in results1:
            if len(set(seq.seq))>1:
                list1.append([seq.id,list(seq.seq)])
            if seq.id=='dm6.' + chrn:
                ref1 = list(seq.seq)
        for seq in results2:
            if len(set(seq.seq))>1:
                list2.append([seq.id,list(seq.seq)])
            if seq.id=='dm6.' + chrn:
                ref2 = list(seq.seq)

        df1 = pd.DataFrame(list1, columns = ['species_chr','seq1'])
        df2 = pd.DataFrame(list2, columns = ['species_chr','seq2'])
        df = pd.merge(df1, df2, how = "inner", on = "species_chr")
        
        #delete chr name from genome.chr name for compatibility with tree
        df.species_chr = df.species_chr.apply(lambda x: x.split(".")[0])
        
        header = pccrid+"\t"+chrn+":"+str(start1)+"-"+str(end1)+"_"+str(start2)+"-"+str(end2)+"_"+strand + "\t" + duplex_header.rstrip("\n")
        if len(df) < 3:
            print("Less than 3 genomes in alignment for PCCR: " + header + " -- skipping")
            continue
        df['seq'] = df.apply(lambda x: x['seq1'] + ["a"]*10 + x['seq2'],1)
        df.drop(['seq1','seq2'], axis = 1, inplace = True)
        #add gaps from reference sequence to structure
        refseq = ref1 + ["."]*10 + ref2
        duplex = list(duplex_left) + ["."]*10 + list(duplex_right)
        try:
            i = 0
            refstruct = []
            for s in refseq:
                if s!= "-":
                    refstruct.append(duplex[i])
                    i+=1
                else:
                    refstruct.append("-")
        except Exception:
            print("ERROR")
            print(header)
            print(refseq)
            print(duplex)
            print("".join(refseq))
            print("".join(duplex))
            print(str(first) + "a"*10 + str(second))
            print("echo -e '"+duplexinp+"' | /home/tools/ViennaRNA/bin/RNAduplex")
            break
#join alignments of panhandles by inner join => sometimes alignment contains gap-only blocks =>
#need to cut gap-only blocks from alignment and its duplex annotation
        arr = np.array(df['seq'].apply(list).tolist())
        arr = np.r_[arr,[refstruct]]
        mask = (arr == "-")
        ix = mask.all(axis=0)
        nogap = np.apply_along_axis("".join, 1, arr[:, ~ix])
        aln = nogap[:-1]
        ann = nogap[-1].replace("(","<").replace(")",">").replace("-",".")
#write stockholm alignments with rnaduplex annotation and id = pccr id + dG of duplex + coordinates  to out.sto file
        out.write("# STOCKHOLM 1.0\n")
        out.write("#=GF AC    " + "".join(pccrid.split(",")[0].split("=")) + "\n")
        out.write("#=GF ID    " + chrn+":"+str(start1)+"-"+str(end1)+"_"+str(start2)+"-"+str(end2)+"_"+strand + "\n")
        out.write("#=GF DE    " + header + "\n")  #header = pccrid+"\t"+chrn+":"+str(start1)+"-"+str(end1)+"_"+str(start2)+"-"+str(end2)+"_"+strand + "\t" + duplex_header.rstrip("\n")
        for seqname, seq in zip (df['species_chr'],aln):
            out.write(seqname + "\t" + seq + "\n")
        out.write("#=GC SS_cons\t" + ann + "\n")
        out.write("//\n")
#write stockholm alignment to temp input file for r-scape
        with open(delete_file + ".sto", "w") as temp_sto:
            temp_sto.write("# STOCKHOLM 1.0\n")
            temp_sto.write("#=GF AC    " + "".join(pccrid.split(",")[0].split("=")) + "\n")
            temp_sto.write("#=GF ID    " + chrn+":"+str(start1)+"-"+str(end1)+"_"+str(start2)+"-"+str(end2)+"_"+strand + "\n")
            temp_sto.write("#=GF DE    " + header + "\n")  #header = pccrid+"\t"+chrn+":"+str(start1)+"-"+str(end1)+"_"+str(start2)+"-"+str(end2)+"_"+strand + "\t" + duplex_header.rstrip("\n")
            for seqname, seq in zip (df['species_chr'],aln):
                temp_sto.write(seqname + "\t" + seq + "\n")
            temp_sto.write("#=GC SS_cons\t" + ann + "\n")
            temp_sto.write("//\n")
#extract subtree from original tree for big maf and write to temp input file for r-scape
        temp_tree = tree.copy("newick")
        temp_tree.prune(list(df.species_chr), preserve_branch_length=True)
        temp_tree.write(format=5, outfile= delete_file + ".nh")
#run R-scape with *.sto and *.nh temp files, write output files to outdir     
#don`t forget <mkdir ./temp> or change it
        process = subprocess.Popen("~/Programs/rscape_v2.0.0.j/bin/R-scape -E 1 -s --nofigures --samplewc --treefile "+ delete_file +".nh --outname " + delete_file + " " + delete_file +".sto",
                                 shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        stdout, stderr = process.communicate()
        out_rscape.write(header + "\n")
        out_rscape.write(stdout.decode())
        if stderr:
            print("ERROR in R-scape " + pccrid + chrn+":"+str(start1)+"-"+str(end1)+"_"+str(start2)+"-"+str(end2)+"_"+strand + ": " + stderr.decode())
        if iteration%1000==0:
            end_time = time.time()
            print(f"done {iteration} PCCRS, time for last 1000: {end_time-start_time} seconds")
            start_time = time.time()
end_time_total = time.time()
print("Time elapsed in total: " + str((end_time_total-start_time_total)/60) + " min")
