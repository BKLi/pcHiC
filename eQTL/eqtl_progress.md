### eQTL/interactions

<!-- toc -->

- [Feb 22 -- IN](#Feb-22----IN)
- [Feb 8 start hfb](#Feb-8-start-hfb)
- [Feb 1 3X gene # sampling](#Feb-1-3X-gene-sampling)
- [Jan 29 1X sampling](#Jan-29-1X-sampling)
- [Jan 28 egenes control](#Jan-28-egenes-control)
- [Jan 25 - Jan 27](#Jan-25---Jan-27)
  * [tss number - controlled sets](#tss-number---controlled-sets)
    + [1. cortical](#1-cortical)
    + [2. hippocampus](#2-hippocampus)
- [Jan 21 / Jan 22](#Jan-21-Jan-22)
  * [Full interaction set](#Full-interaction-set)
- [Jan 18/Jan 19/Jan 20](#Jan-18Jan-19Jan-20)
  * [control groups -- Liver](#control-groups----Liver)
  * [Cortical/Cortex](#CorticalCortex)
    + [Control groups](#Control-groups)
  * [Motor/Spinal Cord](#MotorSpinal-Cord)
- [Jan 17](#Jan-17)
  * [hippocampus control group](#hippocampus-control-group)
- [Jan 16](#Jan-16)
  * [Changing approach to the correct one ... :)](#Changing-approach-to-the-correct-one)
  * [1.1 filter out TSS-eQTL pairs with TSS_distance < 10kb](#11-filter-out-TSS-eQTL-pairs-with-TSS_distance-10kb)
  * [1.2 expand TSS/eQTL coordinates to 5kb bins](#12-expand-TSSeQTL-coordinates-to-5kb-bins)
  * [1.3 intersect](#13-intersect)
  * [1.4 merge](#14-merge)
  * [1.5 extract interaction ID](#15-extract-interaction-ID)
  * [1.6 extract eQTL ID](#16-extract-eQTL-ID)
- [Jan 12](#Jan-12)
- [Jan 11](#Jan-11)
- [Jan 10](#Jan-10)
  * [sampling control eQTL pairs](#sampling-control-eQTL-pairs)
- [Jan 7](#Jan-7)
- [Dec 31](#Dec-31)

<!-- tocstop -->

#### Feb 22

**22599 interneurons.bedpe**
**22598 interneurons_lh_interactions**

* 1981 interneurons_intersect_sig_cortex_merged
* 180 interneurons_intersect_sig_cortex_merged.ID
* 5013 interneurons_intersect_3xcontrol_cortex_merged
* 1698 interneurons_intersect_3xcontrol_cortex_merged.ID
* 1460 interneurons_intersect_1xcontrol_cortex_merged
* 515 interneurons_intersect_1xcontrol_cortex_merged.ID
* 1466 interneurons_intersect_egenes_1x_control_cortex_merged
* 525 interneurons_intersect_egenes_1x_control_cortex_merged.ID
* 1341 interneurons_intersect_sig_cortex_merged.pval
* 1032 interneurons_intersect_1xcontrol_cortex_merged.pval
* 3491 interneurons_intersect_3xcontrol_cortex_merged.pval
* 996 interneurons_intersect_egenes_1xcontrol_cortex_merged.pval
---
* 3037 RGC_intersect_sig_cortex_merged
* 370 RGC_intersect_sig_cortex_merged.ID
* 9525 RGC_intersect_control_cortex_merged
* 2816 RGC_intersect_control_cortex_merged.ID
* 2744 RGC_intersect_1xcontrol_cortex_merged
* 887 RGC_intersect_1xcontrol_cortex_merged.ID

---
* 6382 neurons_intersect_control_cortex_merged
* 1984 neurons_intersect_control_cortex_merged.ID
* 2559 neurons_intersect_sig_cortex_merged
* 269 neurons_intersect_sig_cortex_merged.ID
* 1935 neurons_intersect_1xcontrol_cortex_merged
* 702 neurons_intersect_1xcontrol_cortex_merged.ID

---
* 2564 IPCs_intersect_sig_cortex_merged
* 262 IPCs_intersect_sig_cortex_merged.ID
* 6048 IPCs_intersect_control_cortex_merged
* 1909 IPCs_intersect_control_cortex_merged.ID
* 1804 IPCs_intersect_1xcontrol_cortex_merged
* 604 IPCs_intersect_1xcontrol_cortex_merged.ID

#### Feb 8 start hfb
[GitHub - ijuric/MAPS](https://github.com/ijuric/MAPS)

#### Feb 1 3X gene # sampling
``````bash
shuf -n 17586 Brain_Cortex.allpairs_tss_filtered.genelist > Brain_Cortex.tss_dist_num_filtered.3xgene.genelist
shuf -n 9249 Brain_Hippocampus.allpairs.tss_filtered.genelist > Brain_Hippocampus.allpairs.tss_filtered.3xgene.genelist
``````
* 125518007 Brain_Cortex.allpairs.tss_filt.3xgene.txt
* 65944321 Brain_Hippocampus.tss_filt.3xgene.txt

* 591465 Brain_Hippocampus.allpairs.tss_filt.3xgene.3x.sampled
* 1274736 Brain_Cortex.allpairs.tss_filt.3xgene.3x.sampled

 ``````bash
 awk '{print $6}' Brain_Hippocampus.allpairs.tss_filt.3xgene.3x.sampled_merged | sed 1d | sort | uniq | wc -l
55718
``````

#### Jan 29 1X sampling
* 2494 Brain_Spinal.signif.tss_filtered.genelist

#### Jan 28 egenes control
* 42224282 Brain_Cortex.allpairs.egenes.txt
* 38080655 Brain_Cortex.allpairs.egenes.nonsig.txt

* 22271227 Brain_Hippocampus.allpairs.egenes.txt
* 20144126 Brain_Hippocampus.allpairs.egenes.nonsig.txt

#### Jan 25 - Jan 27
##### tss number - controlled sets
###### 1. cortical
* 5862 Brain_Cortex.signif.tss_filtered.genelist
* 23959 Brain_Cortex.allpairs_tss_filtered.genelist
* 23959 Brain_Cortex.allpairs_tss_filtered.nonsig.genelist

``````bash
 shuf -n 5862 Brain_Cortex.allpairs_tss_filtered.genelist | sort > Brain_Cortex.allpairs_tss_filtered.rand.genelist
 ``````
* 170841224 Brain_Cortex.allpairs_tss_filtered.txt
* 5862 Brain_Cortex.allpairs_tss_filtered.rand.genelist

``````bash
HEADER="$(head -n 1 Brain_Cortex.allpairs.txt)"

grep -f Brain_Cortex.allpairs_tss_filtered.rand.genelist Brain_Cortex.allpairs_tss_filtered.txt > Brain_Cortex.allpairs.tss_dist_num_filtered.txt
sed -i "1 i $HEADER" Brain_Cortex.allpairs.tss_dist_num_filtered.txt 

grep -f Brain_Cortex.allpairs.tss_filtered.nonsig.rand.genelist Brain_Cortex.allpairs_tss_filtered.nonsig.txt > Brain_Cortex.allpairs.tss_dist_num_filtered.nonsig.txt
sed -i "1 i $HEADER" Brain_Cortex.allpairs.tss_dist_num_filtered.nonsig.txt
``````

* 42108651 Brain_Cortex.allpairs.tss_dist_num_filtered.txt
* 39052178 Brain_Cortex.allpairs.tss_dist_num_filtered.nonsig.txt

``````bash
awk '{print $1}' Brain_Cortex.allpairs.tss_dist_num_filtered.txt | sort | uniq | wc -l
5862

awk '{print $1}' Brain_Cortex.tss_dist_num_filtered.nonsig.sampled | sed 1d | sort | uniq | wc -l
5856
``````
* 3X sampling
  * 424912 Brain_Cortex.signif.tss_filtered.bed
  * 1274736 Brain_Cortex.tss_dist_num_filtered.sampled 
  * 1274737 Brain_Cortex.tss_dist_num_filtered.nonsig.sampled
``````bash
 awk '{print $4}' Brain_Cortex.allpairs.tss_dist_num_filtered.sampled.bed | sed 1d | sort | uniq | wc -l
5859
``````
* 9317 cortical_sig.intersect.cortex_tss_dist_num_filtered.sampled_merged.ID
* 9401 cortical_sig_intersect_tss_dist_num_nonsig_sampled_cortex_merged.ID
###### 2. hippocampus
* 3083 Brain_Hippocampus.signif.tss_filtered.genelist
``````bash
awk '{print $2}' Brain_Hippocampus.allpairs.tss_filtered.txt | sed 1d | sort | uniq > Brain_Hippocampus.allpairs.tss_filtered.genelist
``````
* 169283269 Brain_Hippocampus.allpairs.tss_filtered.txt
* 21738109 Brain_Hippocampus.allpairs.tss_num_dist_filtered.txt
``````bash
awk '{print $1}' Brain_Hippocampus.allpairs.tss_num_dist_filtered.txt | sed 1d | sort | uniq | wc -l
3083
``````
* 7481 hippo_intersect_tss_num_dist_filtered_rand_merged.ID
* 7376 Brain_Hippocampus.intersect.allpairs.tss_num_dist_filt.nonsig_merged.ID

#### Jan 21 / Jan 22
##### Full interaction set

* 162359 cortical_full_intersect_sig_cortex_merged.ID

#### Jan 18/Jan 19/Jan 20
##### control groups -- Liver
* 323429 Liver.v7.signif_variant_gene_pairs.txt
* 285390 Liver.signif.bed (TSS > 10kb)
* 35209 hippo_intersect_Liver.signif_merged
* 4561 hippo_intersect_Liver.signif_merged.ID
* 21842 hippo_intersect_Liver.signif_merged_eqtl.ID

##### Cortical/Cortex
* 172489883 Brain_Cortex.allpairs.txt
* 478904 Brain_Cortex.v7.signif_variant_gene_pairs.txt
* 424912 Brain_Cortex.signif.tss_filtered.bed
* 73889 cortical_lh_interactions
* 4165 cortical_intersect_sig_cortex_merged.ID
``````bash
awk '{print $6}' cortical_intersect_sig_cortex_merged |sed 1d | sort | uniq | wc -l
19886
``````
424912 Brain_Cortex.signif.tss_filtered.bed

###### Control groups
1. sampling
``````bash
awk '(NR==1) || ($3 > 10000 || $3 < -10000)' Brain_Cortex.allpairs.txt > Brain_Cortex.allpairs_tss_filtered.txt
``````
* 170841224 Brain_Cortex.allpairs_tss_filtered.txt
* 2124560 Brain_Cortex.rand.tss_filtered.sampled (== 424912\*5)
* 159328628 Brain_Cortex.allpairs_tss_filtered.nonsig.txt
* 201527 cortical_intersect_rand_cortex_merged
* 27933 cortical_intersect_rand_cortex_merged.ID

2. blood
* 8568 cortical_sig_intersect_sig_blood_merged.ID


##### Motor/Spinal Cord

* 1372 motor_intersect_sig_spinal_merged.ID
* 66978 motor_lh_interactions

#### Jan 17
* plot p-value distribution
``````python
hippocampus_significant_eQTL = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.sig.bed", sep="\t")

hippo_sig_intersect_hippo_sig_eQTL_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_intersect_self_sig_eQTL_merged.ID", sep="\t", names=["ID"])

hippo_sig_intersect_hippo_sig_eQTL = pd.merge(hippo_sig_intersect_hippo_sig_eQTL_ID, hippocampus_significant_eQTL, on=["ID"], how="inner")
``````
##### hippocampus control group

**1. sampling random eqtl pairs (sig and non-sig combined, tss > 10kb)**
* 985776 Brain_Hippocampus.tss_filtered.sampled (=197155\*5)
* 144055 hippocampus_tss_filtered_sampled_merged
* 35663 hippocampus_intersect_tss_filtered_sampled_merged.ID
* 90882 hippocampus_tss_filtered_sampled_merged_eQTL_ID

``````bash
 awk '{print $6}' hippocampus_tss_filtered_sampled_merged |sed 1d|sort -k1,1n| uniq > hippocampus_tss_filtered_sampled_merged_eQTL_ID
 ``````

``````python
print(hippo_sig_intersect_hippo_rand["score"].mean())

8.825843871799883
``````
**2. sampling non_sig pairs (tss > 10kb)**
* 170913652 Brain_Hippocampus.allpairs
* filter out pair (TSS < 10kb)
``````bash
awk '(NR==1) || ($3 > 10000 || $3 < -10000)' Brain_Hippocampus.allpairs.txt > Brain_Hippocampus.allpairs.tss_filtered.txt
``````
* 169283269 Brain_Hippocampus.allpairs.tss_filtered.txt
* p-value distribution **(both post-filtering)**

`````python
print(hippocampus_significant_eQTL.groupby(pd.cut(hippocampus_significant_eQTL.pval_nominal, pval_bin)).count()["pval_nominal"])
``````
``````
(1.0000000000000001e-40, 1e-30]       522
(1e-30, 1e-20]                      10686
(1e-20, 1e-10]                      42966
(1e-10, 1e-05]                     128647
(1e-05, 0.1]                        14334
(0.1, 0.5]                              0
(0.5, 1.0]                              0
``````
``````python
print(hippocampus_tss_filtered_eQTL_pval.groupby(pd.cut(hippocampus_tss_filtered_eQTL_pval.pval_nominal, pval_bin)).count()["pval_nominal"])
``````
``````
(1.0000000000000001e-40, 1e-30]         522
(1e-30, 1e-20]                        10686
(1e-20, 1e-10]                        42966
(1e-10, 1e-05]                       130481
(1e-05, 0.1]                       18883137
(0.1, 0.5]                         66748253
(0.5, 1.0]                         83467223
``````
* **cutoff = 0.05**
``````bash
awk '($7 > 0.05)' Brain_Hippocampus.allpairs.tss_filtered.txt > Brain_Hippocampus.allpairs.tss_filtered.nonsig.txt
``````
* 159000155 Brain_Hippocampus.allpairs.tss_filtered.nonsig.txt
* 985776 Brain_Hippocampus.tss_filtered.nonsig.sampled.sorted.bed
* 145441 hippocampus_tss_filtered_nonsig_sampled_merged
* 35574 hippocampus_tss_filtered_nonsig_sampled_merged.ID
* 91461 hippocampus_tss_filtered_nonsig_sampled_merged_eqtl.ID

#### Jan 16
##### Changing approach to the correct one ... :)
**Modified eQTL-preprocess.py; filter out TSS_distance<10kb & add eQTL ID**
**Modified eQTL-merge: adding eqtl_id**

**1. hippocampus**
* Brain_Hippocampus.sig.bed   **221877 (all sig eQTL pairs)**

##### 1.1 filter out TSS-eQTL pairs with TSS_distance < 10kb
* 197155 post-filtering

##### 1.2 expand TSS/eQTL coordinates to 5kb bins
``````bash
awk '{print $1"\t"$2-2500"\t"$2+2500"\t"$4"\t"$5"\t"$6}' Brain_Hippocampus.sig.bed | sed 1d > Brain_Hippocampus.sig.lh.bed
awk '{print $1"\t"$3-2500"\t"$3+2500"\t"$4"\t"$5"\t"$6}' Brain_Hippocampus.sig.bed | sed 1d > Brain_Hippocampus.sig.rh.bed
``````

##### 1.3 intersect
``````bash
bedtools intersect -a Brain_Hippocampus.sig.${INPUT}.bed -b hippocampus_${INPUT}_int\
eractions -wa -wb > hippocampus_intersect_self_sig_${INPUT}
``````
##### 1.4 merge
``````bash
python $HOME/scripts/eQTL-merge.py ${INPUT}_intersect_self_sig_lh ${INPUT}_intersect_self_sig_rh ${INPUT}_intersect_self_sig_merged
``````
* 21210 remained after merging

##### 1.5 extract interaction ID
``````bash
awk '{print $11}' hippocampus_intersect_self_sig_merged |sed 1d|sort|uniq > hippocampus_intersect_self_sig_merged.ID
``````
* 2938 remained after deduplication

``````bash
awk '{print $6}' hippocampus_intersect_self_sig_merged | sed 1d | sort | uniq | wc -l
13245
``````

##### 1.6 extract eQTL ID
``````bash
awk '{print $6}' hippocampus_intersect_self_sig_merged |sed 1d|sort -k1,1n|uniq > hippocampus_intersect_self_sig_eQTL_merged.ID
``````
* 197156 Brain_Hippocampus.sig.bed (# after filtering)
* 13245 hippocampus_intersect_self_sig_eQTL_merged.ID

#### Jan 12 

**II.cortical**
* sig pairs: 
``````python
print(cortical_sig_pairs["pval_nominal"].min())
1.13408e-53
print(cortical_sig_pairs["pval_nominal"].max())
0.00014429200000000002
``````
* filter
``````bash
 awk '(NR==1) || ($4 > 0.00014429200000000002)' Brain_Cortex.allpairs.cut.txt | sort -k1,1 -k2,2n > Brain_Cortex.cut.nonsig.txt
 ``````
 
#### Jan 11

* re-binning distance bins so # of interactions in each bin is uniformly distributed
* Actually.. No need to do that... I think it's fine as long as we have the same propotion of # of interations in each bin, which is independent of how I choose to bin them I think???

**I. hippocampus**
1. significant hippocampus interactions/hippocampus-nonsig-sampled eQTL
* hippocampus: merged/lh = 12026352/17593601
  Uh... proportion seems a bit high?
* hippocampus_sig_sample_nonsig
  * 105135 ID
2. hippocampus_sig_all_sig
  * significant hippocampus interactions/hippocampus-sig-eQTL
  * 38800 ID

**Not really surprisingly, \# of interactions overlapping with eQTL pairs is positively correlated with # of total interactions. Should we control this var?**
  

#### Jan 10

``````bash
awk '{print $10}' hippocampus_intersect_sig_lh | sort | uniq > hippocampus_intersect_sig_single_end.ID

awk '{print $10}' hippocampus_intersect_sig_rh | sort | uniq >> hippocampus_intersect_sig_single_end.ID

sort hippocampus_intersect_sig_single_end.ID | uniq > hippocampus_intersect_sig_single_end.uniq.ID
````````
* \# of IDs : 38800(intersect with both ends) --> 55699(at least one end)
* Both mean and median smaller than that of shared interactions

``````bash
awk '{print $9}' hippocampus_intersect_sig_rh > hippocampus_intersect_sig_single_end.score
awk '{print $9}' hippocampus_intersect_sig_rh >> hippocampus_intersect_sig_single_end.score
``````
* mean and median lower even more

##### sampling control eQTL pairs
```python
hippocampus_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.v7.signif_variant_gene_pairs.txt", sep="\t")

hippocampus_sig_pairs["tss_distance"] = hippocampus_sig_pairs["tss_distance"].abs()

distance_bins = list (np.concatenate((np.arange(0,100000,25000), np.arange (100000,500000,100000), np.arange(500000,1000001,500000))))

hippocampus_sig_pairs.groupby( pd.cut (hippocampus_sig_pairs.tss_distance, distance_bins)).count()["tss_distance"]
```````
```
tss_distance
(0, 25000]           26944
(25000, 50000]       16497
(50000, 75000]       11183
(75000, 100000]       9139
(100000, 200000]     18693
(200000, 300000]      9646
(300000, 400000]      6078
(400000, 500000]      4348
(500000, 1000000]     8080
```
``````python
hippocampus_sig_pairs_binned = hippocampus_sig_pairs.groupby(pd.cut(hippocampus_sig_pairs.tss_distance, distance_bins))
# each sub group correspond to a df with certain tss range

hippocampus_sig_pairs_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _,group in hippocampus_sig_pairs_binned]

# number of samplles in each group of full dataset
sample_size = [group.shape[0]*5 for group in hippocampus_sig_pairs_binned_df_list]
# [134720, 82485, 55915, 45695, 93465, 48230, 30390, 21740, 40400]

# sample_size = [group.shape[0]//5 for group in hippocampus_sig_pairs_binned_df_list] 
test = [hippocampus_sig_pairs_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]

test_concat = pd.concat(test).reset_index().drop(["index"], axis=1)
``````

#### Jan 7

* Some TSS-SNP pairs occur more than once with possibily different p-values (or effect size)? NEED TO FIGURE OUT

#### Dec 31

left hand :$1(chr) $2(start) $3(end)  $10(score) $11(interaction ID)
right hand: $5(chr) $6(start) $7(end) $10(score) $11(interaction ID)

```bash
shuf -n 478904 Brain_Cortex.bed > Brain_Cortex.rand.bed

sort -k1,1 -k2,2n -k5,5 -k6,6n cortical_sig_interaction > cortical_sig_interaction.sorted
awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$11}' cortical_sig_interaction.sorted | sed 1d > cortical_lh_interactions
awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$11}' cortical_sig_interaction.sorted | sed 1d > cortical_rh_interactions
```

