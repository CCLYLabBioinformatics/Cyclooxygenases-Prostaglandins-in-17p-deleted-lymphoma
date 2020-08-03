####bam pathway####
pathway=./QL/sortedByCoord.out.bam
#### get your samples' names (DON'T MODIFY)######
cd $pathway
ls *.bam | while read id ;do 
echo $id ;
done

#### paste your samples' names here ######
shAlox15b.1252_4.sortedByCoord.out.bam
shAlox15b.2865_1.sortedByCoord.out.bam
shAlox15b.2865_2.sortedByCoord.out.bam
shAlox15b.2865_4.sortedByCoord.out.bam
shREN_1.sortedByCoord.out.bam
shREN_2.sortedByCoord.out.bam
shREN_3.sortedByCoord.out.bam
shREN_4.sortedByCoord.out.bam

#### open your case information(DON'T MODIFY) ######
vim $pathway/sample_sampletable1.csv
i
#### write down your case information (DON'T SHOW THIS LINE IN DOCUMENT)######
ids,sample,deal,order
1,shAlox15b.2865_1.sortedByCoord.out.bam,shAlox15b.2865-1,15b,1
2,shAlox15b.2865_2.sortedByCoord.out.bam,shAlox15b.2865-2,15b,2
3,shAlox15b.2865_4.sortedByCoord.out.bam,shAlox15b.2865-4,15b,3
4,sshAlox15b.1252_4.sortedByCoord.out.bam,shAlox15b.1252-4,15b,4
5,shREN_1.sortedByCoord.out.bam,shREN-1,control,5
6,shREN_2.sortedByCoord.out.bam,shREN-2,control,6
7,shREN_3.sortedByCoord.out.bam,shREN-3,control,7
8,shREN_4.sortedByCoord.out.bam,shREN-4,control,8

#### save your case information ######
ESC####按一下你的ESC键
:wq####输入这一行的命令




