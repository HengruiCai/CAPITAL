# CAPITAL: Optimal Subgroup Identification via Constrained Policy Tree Search

This repository is the official implementation of *CAPITAL: Optimal Subgroup Identification via Constrained Policy Tree Search*.

## Introduction

We propose an optimal subgroup selection rule (SSR) that maximizes the number of selected patients, and in the meantime, achieves the pre-specified clinically meaningful mean outcome, such as the average treatment effect. We derive two equivalent theoretical forms of the optimal SSR based on the contrast function that describes the treatment-covariates interaction in the outcome. We further propose a ConstrAined PolIcy Tree seArch aLgorithm (CAPITAL) to find the optimal SSR within the interpretable decision tree class. The proposed method is flexible to handle multiple constraints that penalize the inclusion of patients with negative treatment effect, and to address time to event data using the restricted mean survival time as the clinically interesting mean outcome.  

## Requirements

 - R 3.6
 - `foreach`
 - `doParallel`
 - `policytree`
 - `randomForestSRC`
 - `survival`

## Contents

  1. `README.txt`: implementation details of source code and contents 

  2. Source codes of CAPITAL and data generation environment

     a). `CAPITAL.R`: main function for CAPITAL for Average Treatment Effect with or without Multiple Constraints;

     b). `CAPITAL_surv.R`: main function for CAPITAL for Survival Data with or without Multiple Constraints;
     
     c). `Ture_Calculation.R`: Calculate the true functions and values for Section 5.
     
 

