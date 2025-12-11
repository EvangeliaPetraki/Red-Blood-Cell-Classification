# Red Blood Cell Detection and Fuzzy Classification

This repository contains two MATLAB pipelines for detecting, segmenting, and classifying red blood cells (RBCs) in microscopy images.  
The project implements both traditional image processing and fuzzy logic–based classification.

---

## Contents

- `Prueba_Eva_fuzzy.txt` – Complete RBC detection and **fuzzy logic classifier** implementation.
- `Prueba_Eva_3.txt` – Extended RBC segmentation and refinement workflow (noise removal, reshaping, morphological filtering, rescue logic).

---

## Overview

The project performs:

### **1. Image Preprocessing**
- Reads microscopic RBC images.
- Works primarily on the **green channel**, extracted due to best contrast for RBC boundaries.
- Computes the histogram and applies a custom **thresholding method** based on the first peak half-maximum (FWHM) to generate an initial binary mask.

### **2. Segmentation and Morphological Processing**
Both scripts apply:
- Noise removal through erosion and dilation.
- Region labeling (`bwlabel`) to count and track individual cells.
- Removal of small objects based on adaptive thresholds.
- Extraction of geometric descriptors:
  - **Area**
  - **Eccentricity**
  - **Solidity**
  - **Bounding boxes**

The second script includes additional:
- Border-touching removal
- Resizing/rescue pipeline for borderline or ambiguous cells
- Solidity-based re-evaluation after erosion
- Final reconstruction of a cleaned binary mask

### **3. Feature-Based Classification**
The fuzzy classifier (from `Prueba_Eva_fuzzy.txt`) defines:

- **Inputs:**
  - Area  
  - Solidity  
  - Eccentricity  

- **Output:**  
  - 0 → None  
  - 1 → Normal RBC  
  - 2 → Pathological RBC  

Membership functions are defined using trapezoidal (`trapmf`) and triangular (`trimf`) shapes.  
Classification is done through a rule base (17+ rules) combining feature conditions.

### **4. Visual Annotation**
The pipeline overlays bounding boxes onto the original image:

- **Red boxes:** normal RBCs  
- **Green boxes:** pathological RBCs  

It prints the total number of cells and the number of pathological detections.

---

## Usage

1. Place your image files in the same directory as the scripts.
2. Update the `imread(...)` lines to load the desired microscopy image.
3. Run the `.m` script in MATLAB.

Both scripts produce:
- Intermediate segmentation visualizations
- Plots of area, eccentricity, and solidity distributions
- A final annotated image showing all detected and classified cells

---

## Requirements

- MATLAB R2018a or newer
- Image Processing Toolbox
- Fuzzy Logic Toolbox

---

## Citation

This work was accepted and presented as a full paper at the International Joint Conference on Computational Intelligence 2025.

If you use this code, please cite the corresponding publication (details to be added).


