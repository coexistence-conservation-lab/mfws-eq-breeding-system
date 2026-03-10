# Title

Here we present the datasets and workflows used for Brockett *et al*., (in review) *Factors affecting multiple paternity: insights from the eastern quoll (*Dasyurus viverrinus*)*. Related datasets are also available in the [ANU data repository](https://10.25911/50fr-x765).

## Highlights 

  - Eastern quolls exhibited multiple paternity in ≥47% of litters
  - Reproductive success was greater for younger, heavier, and dark-coloured individuals 
  - High rates of multiple paternity in Dasyuridae were correlated with shorter male lifespan, but not female lifespan, sexual size-dimorphism, or aggression

## Abstract

Breeding systems influence population dynamics, effective population size, and genetic diversity, with particularly strong consequences for small, isolated, or reintroduced populations. To examine both population-level drivers and broader evolutionary constraints on multiple paternity, **we combined long-term demographic data and 1,745 single nucleotide polymorphisms from a reintroduced population of eastern quolls** (*Dasyurus viverrinus*) with a comparative analysis across the family Dasyuridae. 

Pedigree reconstruction revealed **high rates of multiple paternity**, with 47–85% of litters sired by more than one male depending on litter size, and up to three sires per litter. **Individual reproductive success was significantly associated with body weight, age, and colour morph**, indicating that phenotypic variation can lead to unequal genetic contributions within populations. Across the Dasyuridae family, shorter lifespans were associated with higher rates of multiple paternity, while sexual size-dimorphism and intersex aggression were not.

These findings indicate that while individual traits shape reproductive outcomes within populations, **the prevalence of  multiple paternity is more strongly influenced by species-level life-history strategies**. Incorporating mating system dynamics into conservation planning may therefore refine predictions of genetic outcomes in reintroduced or conservation-dependent species.

## Approach

We posed the question: *"What affects breeding success in  population"* in a hurdle-model approach.
  
  1. Did **morphological** (e.g., pelage colour morph) or **genetic** (e.g., heterozygosity) factors drive breeding success? 
  2. For individuals that breed successfully, do these factors influence the number of offspring produced?

Note that in this analysis "successful" reproduction is equivalent to recruiting offspring into the next generation, rather than conception/birth being the endpoint. This is because individuals can only be captured and  microchipped/sampled if they survive from birth to dispersal from the natal den (and, therefore, they are captured in population monitoring).

## Repository structure

This repository follows an organised structure for clarity and reproducibility:

### Files

  - `.gitattributes`: repository-specific Git settings 
  - `.gitignore`: lists files and filetypes to be excluded from version control 
  - `LICENSE`: MIT license 
  - `README.md`: project overview 
  - `analyses.Rmd`: R Markdown with analysis workflows 
  - `functions.R`: project-specific functions sourced by `analyses.Rmd`
  - `project.Rproj`: RStudio project file for consistent setup 
  - `tutorial.html`: rendered summary of analyses and results 

### Folders

  - `archive\`: superseded datasets and workflows (if applicable)
  - `input\`: raw and reference datasets (e.g., geospatial files, spreadsheets) 
  - `output\`: output from analyses (e.g., maps, plots, spreadsheets) 
  - `metadata\`: API keys and citation files 

## Licence

Unless otherwise stated, all code in this repository is licensed under the MIT License. We kindly ask that you cite the relevant publication(s) and this repos if you reuse or adapt our code. 