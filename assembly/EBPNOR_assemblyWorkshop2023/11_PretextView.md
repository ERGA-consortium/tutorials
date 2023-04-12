# PretextView tutorial

There are many ways to view a Hi-C contact map, but today we are going to use PretextView. To download the PretextView desktop application, click [here](https://github.com/wtsi-hpag/PretextView/releases), and pick a release that is suitable for your laptop. 

## How to curate your assemblies

### Step 1: Tweaking your settings

 To start the curation process, open the PretextView desktop application, click **load map**, and navigate to where you saved your pretext files. 
 
 When you load your Hi-C contact map, this is what it should look like, depending on which colour scheme you have chosen. For this tutorial we will be using **"Blue-Orange Divergent"**.  
 
 <img width="962" alt="Screenshot 2023-02-06 at 13 43 38" src="https://user-images.githubusercontent.com/110542053/216974988-510ca53d-a3f9-4d84-8a0f-77ce8d53057a.png">

 
 Here, each square represents a scaffold (which after curation will hopefully all be of chromosome length). The red diagonal line shows where the strongest contact signals are between the DNA sequences. 

 Go ahead and look at the extensions for your contact maps. These overlays makes it easier to figure out where the unplaced and wrongly oriented scaffolds are supposed to go. 

 For this tutorial, turn on the **Gaps** extension, and turn the **Gamma Min** and **Gamma Mid** sliders down to zero, and the **Gamma Max** slider all the way up. This will make it easier to see where there are gaps in the assembly, and increase the contrast so the contact signals will be easier to interpret. 

 ### Step 2: Moving scaffolds in PretextView

For this step you´ll need a computer mouse. Before you do your first edit, try to orient yourself in the Hi-C contact map. Pressing **U** removes the main menu. Don´t worry, pressing **U** again brings it right back if you want to change any settings. 

Try to move around in the contact map. Scrolling the mouse wheel zooms you in and out. You´ll notice that the zoom is not affected by where your cursor is, so how do we zoom in on one particular scaffold? By dragging the map and placing the area you want to look at in the centre of the screen. To drag the map, click and hold the right mouse button, and move the mouse around. Move around the map for a bit, and look at the scaffolds. Do some look fragmented? Move to the far right bottom corner of the contact map. Are there many smaller scaffolds there?

When your comfortable with this movement, let´s try to bring up some of the other menus! Pressing **E** activates "Edit mode". When this mode is active, the cursor changes. Do you see that the further away from the red diagonal you move, the larger an area of the scaffold is marked by green? This green indicator shows which part of the scaffold you are picking up. Try to make an edit! Cut out a chunk of a scaffold and move it someplace else in the contact map. Don´t worry about whether the edit is correct, we´ll delete all the test edits before we start the proper curation process. 

Move to the far right of the contact map. Are there any smaller, unplaced scaffolds with clear, red contact signals that you think would go well in any of the larger scaffolds? Hover over them, press the space bar, and move the cursor without clicking the left button. When you have oriented yourself to where you want to place the unplaced scaffold, click the left cursor. If the piece fits best on the end of one of the larger squares, press **S** while in Edit mode to toggle the Snap function. Do you notice that the unplaced scaffold "travels" differently across the contact map, in a skipping motion? This lets you "snap" the unplaced scaffold in place at the end of the larger scaffolds. 

Now that you know how to move around and make edits in PretextView, delete all your edits with **Q** while in Edit mode, and press **U** to check to see that all your test edits are gone. You are now ready to edit the assembly for real!

### Step 3: Editing your TPF

When editing a Hi-C contact map, you need to make *the same* edits to the TPF-file you generated with the Rapid curation suite. 

To make edits to your TPF in the command line, type:

```
nano filename.tpf
```

Here are some different scenarios that you may encounter:

#### "I want to add an unplaced scaffold to the end of one of the larger, chromosome sized scaffolds"

No worries, just enter Edit mode, toggle the snap function, and snap the unplaced scaffold onto the end of the larger chromosome sized scaffold. No need to edit the TPF!


#### "I want to add an unplaced scaffold to the middle of one of the larger, chromosome sized scaffolds, where there is a gap"

This can be done! Just grab the unplaced scaffold with the space bar and move it to the gap (this is where the "Gaps" extension comes in handy). Press **E** to exit edit mode, and press **U** to see your edits. Here you can see the exact position you moved the scaffold to. 

Then, go to the TPF, and add a **>** to the line you want to remove. For instance, say that we have a gap in scaffold 1, at position 408378. Before editing, the TPF will look like this:

```
?	scaffold_1:1-408378     scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:408579-1794055	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:1794256-2710345	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:2710546-2734866	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:2735067-3369476	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:3369677-3660413	scaffold_1	PLUS
```

After editing, your TPF should look like this:

```
?	scaffold_1:1-408378     scaffold_1	PLUS
>GAP     TYPE-2  200
?	scaffold_1:408579-1794055	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:1794256-2710345	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:2710546-2734866	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:2735067-3369476	scaffold_1	PLUS
GAP     TYPE-2  200
?	scaffold_1:3369677-3660413	scaffold_1	PLUS
```

You have made your first edit! There is no need to do anything to the scaffold you have now placed. 


#### "I want to add an unplaced scaffold to the middle of one of the larger, chromosome sized scaffolds, where there is NOT a gap"

While this is possible, we will not be doing those kind of edits in this course. You need a higher resolution editor to be able to create gaps in the correct place in the assembly (such as HiGlass, which you can read about [here](http://higlass.io/)).

#### "Some of my unplaced scaffolds have ambiguous contact signals, what do I do?"

Ask us for help! You can also read more about these ambiguous signals in GRIT´s documentation [here](https://gitlab.com/wtsi-grit/rapid-curation/-/blob/main/PretextView%20-%20Tutorial.pdf). 

#### "I want to make edits within one of the larger scaffolds, what do I do"

In some instances you´ll want to invert segments in your chromosome sized scaffolds. If there are gaps in the breakpoints where you want to make the inversion, you simply edit the TPF as shown above. However, if there are no gaps, do not make the edit. To invert segments, press spacebar again after you have picked up the segment you want to invert. 

### Step 4: Painting your scaffolds

You have made your edits, and now you hopefully have eight large scaffolds matching *Athalia rosae´s* karyotype. Before finishing your assembly, you need to "paint" them. What does this mean? You need to mark which scaffolds are part of the same "super-scaffolds" or chromosomes, so they´ll all have the same name in the final FASTA. 

Enter the "Scaffold Edit Mode" by pressing **S**. Go to the bottom right corner, and left click on the smallest scaffold you want to include as a chromosome. 

While clicking, hold **A**, and drag in a diagonal line (following the contact signal) till you reach the end of the chromosome. Repeat until you have reached the end in the left top corner, and have painted all the scaffolds. 


### Step 5: Finishing your assembly

To finish your assembly you need to:

#### 1. Create an AGP file

Press **U** to bring up the main menu. Press the "Generate AGP" button, and create a out.pretext.agp file. 

#### 2. Transfer your out.pretext.agp file to saga

Bring the AGP file back to your saga working directory by opening another terminal window, logging into saga. and going to the folder where you saved your TPF-file. When you are in the right location, use the code below to copy your file:

```
scp -r out.pretext.agp <username>@saga.sigma2.no:/cluster/projects/nn9984k/folder_with_tpf
```

Enter your password and the file will be transferred to the directory. 

#### 3. Create a new TPF and new FASTA from the original edited TPF and AGP file

Activate the rapid_curation conda environment. Using the rapid_pretext2tpf_XL.py script, combine the edited TPF and the generated AGP to create a new TPF. 

```
eval "$(/cluster/projects/nn9984k/miniconda3/bin/conda shell.bash hook)" 

conda activate curation

# create new tpf

python /cluster/projects/nn9984k/opt/rapid-curation/rapid_pretext2tpf_XL.py \
iyAthRosa_clean.fa.tpf \
out.pretext.agp_1
```

Using the new TPF, the chrs.csv-file generated with the Rapid curation suite, and the original fasta, you can create a new fasta with the rapid_join.pl script:


```
# create new fasta

perl /cluster/projects/nn9984k/opt/rapid-curation/rapid_join.pl -fa iyAthRosa_clean.fa \
-tpf rapid_prtxt_XL.tpf \
-csv chrs.csv \
-out hap1_clean.fasta
```

And now you are left with a complete, curated, haplotype resolved whole-genome assembly! Congratulations, you champ!

![leo_gif](https://user-images.githubusercontent.com/110542053/206199166-141c2f3b-2f9c-42a4-913f-62e3913511fa.gif)


|[Previous](https://github.com/ebp-nor/genome-assembly-workshop-2023/blob/main/10_Rapid_curation.md)|
|---|
