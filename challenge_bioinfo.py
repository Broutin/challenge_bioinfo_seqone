import pysam

#ouvrir fichier
bamfile = pysam.AlignmentFile("wetransfer-ac920e/HD200-10ng-E02A-B04_S16.bam", "rb")
ref=pysam.FastaFile("human_g1k_v37.fasta")

liste_chromosomes = bamfile.references


#Initialisation variable
nb_align=0
nb_match=0
ratio=0


for chromosome in liste_chromosomes :
    #dans le fichier genome_ref M=MT
    if chromosome =="M":
        chromosome_ref="MT"
    else :
        chromosome_ref=chromosome
        
    reference=ref.fetch(chromosome_ref)
    for read in bamfile.fetch(chromosome):
        if not read.is_unmapped: #tous ceux qui sont alignes

            #caculer Nombre de bp alignees et identique a la reference
            tuple_read=read.get_aligned_pairs(with_seq=True)
            for tuple in tuple_read :
                if tuple[1] !=None: # Ne prend que ceux alignees sur la ref
                    if reference[tuple[1]]==tuple[2]:
                        nb_match+=1
            
            #calculer Nombre de bp alignees sur la reference (M)
            for cigar in read.cigartuples :
                if cigar[0] == 0 : 
                    nb_align+=cigar[1]
print(nb_align)
print(nb_match)


#ratio
ratio=float(nb_match)/nb_align
print(ratio)

#fermer fichiers
bamfile.close()
ref.close()