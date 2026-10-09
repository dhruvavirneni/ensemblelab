# Analysis
* incomplete documentation section...just contains rough notes

## Central question: How does increasing conformer sampling affect the structural and energetic characterization of a molecule, and how can those changes help determine optimized conformer sampling sizes in ensemble analysis?
 * other stuff i'm wondering: in a given Ensemble with n conformer sample size, will all metrics converge relatively at the same n?

## goals and stuff
Investigate how many conformers are needed to adequately characterize a molecule's conformational ensemble, and whether that number can be predicted from the molecule's properties. The first stage is to study how ensemble-level metrics change as the number of generated conformers increases. Starting with a set of molecules spanning different sizes and levels of structural complexity, I will generate conformational ensembles at progressively larger sampling sizes, optimize their geometries, and standardize the resulting structures and energies for consistent comparison.

Ensemble wide metrics worth calculating/investing convergence on:
 * relative energy distribution: (for each conformer, calculate energy relative to lowest energy conformer) --> good for comparing conformational stability
 * boltzmann populations: assign each conformer a statistical weight based on its energy and a chosen temperature, then normalize the weights into estimated populations --> tells us which conformers contribute most to population weighted ensemble
 * conformational entropy and effective ensemble size: calculate Shannon entropy to measure how 'spread out' the ensemble is --> whether an ensemble is dominated by a few conformers or distribute across many...remember to calculate populations at cluster level too
    * according to google: calculating populations at the conformer cluster level is necessary because Shannon entropy is highly sensitive to how states are divided. If your ensemble contains many near-duplicate structures, treating them as distinct conformers artificially inflates the entropy $H$ and distorts the effective ensemble size ($N_{eff}$)
 * torsional state diversity: track distribution of torsions around rotatable bonds --> which internal rotations account for conformational variations
  * can track by simply tracking dihedral angles around given bonds, or assigning torsions to their 'discrete states'(gauche+, gauche-, anti) **investigate if one will be easier for later ML use cases**
 * radius of gyration
  * how spread out atoms are around conformer's center of mass --> whether ensemble is more compact or extended --> output mean, std, range, etc.
 * ensemble molecular properties: solvent accessible surface area, dipole moment, etc...intramolecular hbonds --> whether conformational changes affect properties beyond solely geometry

 * another cool thing (suggested by AI): cluster discovery rate: num of newly discovered clusters per additional batch of conformers --> could complement minimum energy and structural convergence curves


For each sampling size, I will calculate metrics such as minimum energy, pairwise RMSD, structural clustering, and other measures of conformational diversity. Plotting these metrics against the number of conformers should reveal how quickly each stabilizes and whether additional sampling continues to uncover meaningful structural or energetic differences. Since different metrics may converge at different rates, I will investigate how to define both metric-specific convergence points and an overall sampling size that satisfies the chosen criteria.


## pipeline + plan
 1. define experiement 
    * pick molecule of choice, sampling size tests, random seeds, optimization settings, etc.
 2. generate conformers
    * generate large reference pool in Ensemble...use same random seed to generate ensembles with different nconformers
 3. optimize and standardize conformers
    * consistent optimization, deduplication, etc.
 4. analyze each sampling size
 5. quantify convergence per sample size
    * measure metric changes between sample sizes (metrics to look at listed above)
 6. train ML model with generated dataset of multiple molecules, having run this analysis
 7. test ML model, using molecule input and other quick molecular metrics like logP, etc. etc. as parameters and outputting optimal n (sample size) for effectively using this molecule for ensemble analysis
 8. maybe build a web app interface or something that maintains records of ML analyses conducted and runs new ones to output n if an inputted molecule hasn't been analyzed by anyone yet (much later)
