# GWAC

GWAC is a Google Earth Engine (GEE) module for atmospheric correction of water targets from satellite imagery.

This repository currently includes implementations for:

- Landsat 8
- Landsat 9
- Sentinel-2

Rayleigh scattering correction is performed with pre-trained **GRAYCO** models loaded from Earth Engine assets.

## Overview

GWAC is designed for water-target atmospheric correction workflows over inland and coastal waters.  
The current implementations combine:

- ancillary atmospheric data attachment
- water masking and QA filtering
- pressure-scaled Rayleigh optical thickness
- absorbing gas (ozone) correction
- GRAYCO-based Rayleigh correction
- SWIR-based aerosol extrapolation
- remote sensing reflectance (Rrs) retrieval

## Repository structure

```text
GWAC/
├── GWAC_L8_module.js
├── GWAC_L9_module.js
├── GWAC_S2_module.js
├── GWAC_example.js
├── README.md
└── LICENSE
```

## Contact
If you use this code in research, please cite the related paper or project documentation once available.

## Contact
Questions, suggestions, and feedback are always welcome! Please feel free to open an issue, or contact me by email: dianazhao0105@gmail.com
