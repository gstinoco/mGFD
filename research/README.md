# mGFD Research Laboratory :microscope:

<div align="center">

<img src="docs/logo/logo.png" alt="mGFD logo" width="680" style="margin: 20px 0;">

**Experimental tools, benchmarks, visualizations, and reference datasets for the mGFD method.**

### :link: Quick Links
[![🧮 Methodology](https://img.shields.io/badge/🧮-Methodology-purple)](METHODOLOGY.md) [![🗂️ Datasets](https://img.shields.io/badge/🗂️-Datasets-blue)](DATASETS.md) [![🎬 Results](https://img.shields.io/badge/🎬-Results-purple)](RESULTS.md) [![👥 Team](https://img.shields.io/badge/👥-Research%20Team-blue)](#scientist-research-team) [![🏭 Partners](https://img.shields.io/badge/🏭-Industry%20Partners-0B1B3A)](#factory-industry-partners--sponsorships)

</div>

---

## :clipboard: Table of Contents
- [📚 Scientific Documentation](#books-scientific-documentation)
- [📈 Experiment Reproduction](#chart_with_upwards_trend-experiment-reproduction)
- [📊 Computational Performance Profile](#bar_chart-computational-performance-profile)
- [👥 Research Team](#scientist-research-team)
- [🏭 Industry Partners & Sponsorships](#factory-industry-partners--sponsorships)
- [📄 Scientific References](#books-scientific-references)

---

This directory contains the entire suite of experimental tools, benchmarks, visualizations, and original datasets (real lakes) that support the scientific publication of the mGFD method.

Unlike the project root (which functions as the installable Python package), this folder (`/research/`) is the **research laboratory** where scripts consume the `mGFD` library to generate geographic studies and numerical stability analyses.

---

## :books: Scientific Documentation

To maintain readability, the deep scientific details of the laboratory are divided into specialized markdown files. Please visit the following pages for comprehensive information:

*   **[🧮 Mathematical Methodology (METHODOLOGY.md)](METHODOLOGY.md)**: Detailed derivations of the mGFD scheme, Taylor expansions, consistency conditions, and the Upwind stabilization scheme.
*   **[🗂️ Benchmark Datasets (DATASETS.md)](DATASETS.md)**: A complete catalog of the 20 real-world geographic point clouds, density scaling (1-5), and KD-Tree stencil structures.
*   **[🎬 Experimental Results (RESULTS.md)](RESULTS.md)**: Visual demonstrations (GIFs/PNGs) and error bounds for the Poisson, Heat, Wave, and Advection-Diffusion equations.

---

## :chart_with_upwards_trend: Experiment Reproduction

To reproduce the results generated in the scientific publications, or to generate massive performance benchmarks across all available point clouds, you can use the centralized parameter sweep orchestrator:

```bash
# Move to the research/codes directory
cd research/codes

# Run all experiment batches at once
python sweep.py
```

### ⚙️ High-Performance Configuration (`sweep_config.json`)

The orchestrator reads `codes/sweep_config.json` to generate a combinatorial execution matrix. You can modify this file to benchmark the new **GPU Acceleration** and **Matrix-Free** capabilities introduced in v0.11.0:

```json
{
    "nvec": [12, 16, 20],
    "linear_solver": ["spsolve", "bicgstab", "gmres"],
    "device": ["cpu", "cuda"],
    "matrix_free": [false, true],
    "preconditioner": [null, "ilu"],
    "save": true,
    "upwind": [true, false]
}
```

> [!TIP]
> **Benchmarking GPU vs CPU:** Simply add `"cuda"` to the `"device"` array. The sweep will automatically run the problem on the CPU, and then re-run it on the GPU, outputting CSVs that allow for direct speedup comparisons.

The raw metrics are exported as JSON files inside `results/`. **However, for your convenience, `sweep.py` will automatically crawl all generated JSONs at the end of the execution and compile a master `sweep_summary_YYYYMMDD_HHMMSS.csv` inside the `results/` folder.** This file contains a tabular view of all combinations (Equation, Scale, Device, Matrix-Free, Time, Error, etc.), making it incredibly easy to plot and analyze CPU vs GPU performance.

---

## :bar_chart: Computational Performance Profile

The `mGFD` suite is rigorously benchmarked. The following visual profile and data table summarize **500 automated tests** evaluating 5 families of equations across 20 distinct geographic point clouds and 5 density scales.

<div align="center">
  <img src="docs/images/benchmark_profile.png" alt="Benchmark Profile" width="800"><br/>
  <sub>Averaged execution time and memory footprint for standard unstructured domains</sub>
</div>

### :bar_chart: Global Averages (Summary)

| Equation | :stopwatch: Avg. Time (s) | :floppy_disk: Avg. Memory (MB) | :triangular_ruler: Avg. Error (RMSE) |
|:---|:---:|:---:|:---:|
| **Poisson** | `0.230` | `3.2` | `1.86e-06` |
| **Advection** | `0.934` | `140.6` | `7.63e-05` |
| **Advection-Diffusion** | `0.923` | `127.0` | `8.22e-04` |
| **Heat** | `2.425` | `137.0` | `1.17e-07` |
| **Wave** | `2.434` | `112.4` | `3.18e-06` |

*Note: Computations rely on sparse KD-Tree stencils and the `scipy.sparse.linalg.bicgstab` solver.*

## :scientist: Research Team

<div align="center">

### :star2: Meet the Team
*Researchers and graduate students advancing meshless computational methods*

</div>

### :busts_in_silhouette: Main Researchers

<table align="center" width="100%" cellspacing="14">
  <tr>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/gtinoco.webp" width="100" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Dr. Gerardo Tinoco Guerrero</b> :mexico:<br/>
        <sub>Numerical Methods & Computational Mathematics</sub><br/><br/>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:gerardo.tinoco@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
        <a href="https://orcid.org/0000-0003-3119-770X"><img alt="ORCID" src="https://img.shields.io/badge/ORCID-0000--0003--3119--770X-green"></a>
        <a href="https://www.researchgate.net/profile/Gerardo-Tinoco-Guerrero"><img alt="ResearchGate" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </div>
    </td>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/dmota.webp" width="100" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Dr. Francisco Javier Domínguez Mota</b> :mexico:<br/>
        <sub>Applied Mathematics & Finite Difference Methods</sub><br/><br/>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:francisco.mota@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
        <a href="https://orcid.org/0000-0001-6837-172X"><img alt="ORCID" src="https://img.shields.io/badge/ORCID-0000--0001--6837--172X-green"></a>
        <a href="https://www.researchgate.net/profile/Francisco-Dominguez-Mota"><img alt="ResearchGate" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </div>
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/jagt.webp" width="100" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Dr. José Alberto Guzmán Torres</b> :mexico:<br/>
        <sub>Engineering Applications & Artificial Intelligence</sub><br/><br/>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:jose.alberto.guzman@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
        <a href="https://orcid.org/0000-0002-9309-9390"><img alt="ORCID" src="https://img.shields.io/badge/ORCID-0000--0002--9309--9390-green"></a>
        <a href="https://www.researchgate.net/profile/Jose-Guzman-Torres"><img alt="ResearchGate" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </div>
    </td>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/harias.webp" width="100" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Dr. Heriberto Árias Rojas</b> :mexico:<br/>
        <sub>Engineering Applications</sub><br/><br/>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:heriberto.arias@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
        <a href="https://orcid.org/0000-0002-7641-8310"><img alt="ORCID" src="https://img.shields.io/badge/ORCID-0000--0002--7641--8310-green"></a>
        <a href="https://www.researchgate.net/profile/Heriberto-Arias-Rojas"><img alt="ResearchGate" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </div>
    </td>
  </tr>
</table>

### :mortar_board: Ph.D. Research Students

<table align="center" width="100%" cellspacing="14">
  <tr>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/gpj.webp" width="80" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Gabriela Pedraza-Jiménez</b><br/>
        <img alt="Ph.D. Research Student" src="https://img.shields.io/badge/Ph.D.-Research%20Student-2E8B57?style=flat-square"><br/><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:2220157h@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </div>
    </td>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/eci.webp" width="80" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Eli Chagolla-Inzunza</b><br/>
        <img alt="Ph.D. Research Student" src="https://img.shields.io/badge/Ph.D.-Research%20Student-2E8B57?style=flat-square"><br/><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:1137626b@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </div>
    </td>
  </tr>
</table>

### :mortar_board: M.Sc. Research Students

<table align="center" width="100%" cellspacing="14">
  <tr>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/jlgf.webp" width="80" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Jorge L. González-Figueroa</b><br/>
        <img alt="M.Sc. Research Student" src="https://img.shields.io/badge/M.Sc.-Research%20Student-green?style=flat-square"><br/><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:1718717h@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </div>
    </td>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/cnmb.webp" width="80" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Christopher N. Magaña-Barocio</b><br/>
        <img alt="M.Sc. Research Student" src="https://img.shields.io/badge/M.Sc.-Research%20Student-green?style=flat-square"><br/><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:1339846k@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </div>
    </td>
  </tr>
</table>

### :mortar_board: Undergraduate Research Students

<table align="center" width="50%" cellspacing="14">
  <tr>
    <td width="100%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/team/mgfl.webp" width="80" style="border-radius: 50%; margin-bottom: 15px;"><br/>
        <b>Maria Goretti Fraga-Lopez</b><br/>
        <img alt="Undergraduate Research Student" src="https://img.shields.io/badge/Undergrad-Research%20Student-green?style=flat-square"><br/><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a><br/><br/>
        <a href="mailto:1702174b@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </div>
    </td>
  </tr>
</table>

</div>

---

## :factory: Industry Partners & Sponsorships

<div align="center">

*Collaboration between academia, government, and industry to accelerate real-world impact*

<table align="center" width="100%" cellspacing="14">
  <tr>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/partners/secihti.webp" alt="SECIHTI Logo" height="70" style="margin-bottom: 15px;"><br/>
        <b>🔬 SECIHTI</b><br/>
        <sub>Science, Humanities, Technology & Innovation</sub><br/><br/>
        <a href="https://conahcyt.mx/"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-0B1B3A?style=flat-square"></a>
        <img alt="Type: Government" src="https://img.shields.io/badge/🏛️%20Type-Government-1A3A6B?style=flat-square">
        <br/><br/>
        <b>Role:</b> Financial support and research grants<br/>(Project 1022677)
      </div>
    </td>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/partners/siiia.webp" alt="SIIIA MATH Logo" height="70" style="margin-bottom: 15px;"><br/>
        <b>🏭 SIIIA MATH</b><br/>
        <sub>Soluciones de Ingeniería, México</sub><br/><br/>
        <a href="http://www.siiia.com.mx"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-blue?style=flat-square"></a>
        <img alt="Type: R&D Company" src="https://img.shields.io/badge/📊%20Type-R%26D%20Company-orange?style=flat-square">
        <br/><br/>
        <b>Role:</b> Technology transfer, applied R&D, and computational infrastructure
      </div>
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/partners/umsnh.webp" alt="UMSNH Logo" height="70" style="margin-bottom: 15px;"><br/>
        <b>🎓 UMSNH</b><br/>
        <sub>Universidad Michoacana de S.N.H.</sub><br/><br/>
        <a href="http://www.umich.mx"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-darkred?style=flat-square"></a>
        <img alt="Type: University" src="https://img.shields.io/badge/🏷️%20Type-University-1A3A6B?style=flat-square">
        <br/><br/>
        <b>Role:</b> Academic collaboration, research infrastructure, and graduate formation
      </div>
    </td>
    <td width="50%" valign="top" align="center">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 20px; height: 100%;">
        <img src="docs/partners/cimne.webp" alt="CIMNE Logo" height="70" style="margin-bottom: 15px;"><br/>
        <b>🌿 Aula CIMNE-Morelia</b><br/>
        <sub>Research collaboration space</sub><br/><br/>
        <a href="https://aulas.cimne.com/aula/aula-morelia/"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-orange?style=flat-square"></a>
        <img alt="Type: Research Center" src="https://img.shields.io/badge/🧮%20Type-Research%20Center-EE9B00?style=flat-square">
        <br/><br/>
        <b>Role:</b> Numerical methods environment and academic-industry collaboration
      </div>
    </td>
  </tr>
</table>

</div>

---

## :books: Scientific References

If you use this laboratory or the **mGFD method** in your research, please cite the primary mathematical formulation:

1. **Tinoco-Guerrero, G.**, Domínguez-Mota, F. J., Pedraza-Jiménez, G., Guzmán-Torres, J. A., & Tinoco-Ruiz, J. G. (2025). *"mGFD: A meshless generalized finite difference method"*. **Computers and Mathematics with Applications**, 195, 396-418. [DOI: 10.1016/j.camwa.2025.07.034](https://doi.org/10.1016/j.camwa.2025.07.034)

### Related Publications

2. **Pedraza-Jiménez, G.**, Tinoco-Guerrero, G., Domínguez-Mota, F. J., Guzmán-Torres, J. A., & Tinoco-Ruiz, J. G. (2025). *"mGFD: CloudGenerator"*. **Software Impacts**, 23, 100721. [DOI: 10.1016/j.simpa.2024.100721](https://doi.org/10.1016/j.simpa.2024.100721)
3. **Tinoco-Guerrero, G.**, Arias-Rojas, H., Guzmán-Torres, J. A., Román-Gutiérrez, R., & Tinoco-Ruiz, J. G. (2023). *"A meshless finite difference scheme applied to the numerical solution of wave equation in highly irregular space regions"*. **Computers & Mathematics with Applications**, 136, 25-33. [DOI: 10.1016/j.camwa.2023.01.035](https://doi.org/10.1016/j.camwa.2023.01.035)
4. **Tinoco-Guerrero, G.**, Domínguez-Mota, F. J., Guzmán-Torres, J. A., Román-Gutiérrez, R., & Tinoco-Ruiz, J. G. (2023). *"Study of the stability of a meshless generalized finite difference scheme applied to the wave equation"*. **Frontiers in Applied Mathematics and Statistics**, 9, 1214022. [DOI: 10.3389/fams.2023.1214022](https://doi.org/10.3389/fams.2023.1214022)
5. **Tinoco-Guerrero, G.**, Domínguez-Mota, F. J., Guzmán-Torres, J. A., & Tinoco-Ruiz, J. G. (2022). *"Numerical solution of diffusion equation using a method of lines and generalized finite differences"*. **Revista Internacional de Métodos Numéricos para Cálculo y Diseño en Ingeniería**, 38(2), 24. [DOI: 10.23967/j.rimni.2022.06.003](https://doi.org/10.23967/j.rimni.2022.06.003)

### BibTeX Citation

```bibtex
@article{tinoco2025mgfd,
  title={mGFD: A meshless generalized finite difference method},
  author={Tinoco-Guerrero, Gerardo and Domínguez-Mota, Francisco Javier and Pedraza-Jiménez, Gabriela and Guzmán-Torres, José Alberto and Tinoco-Ruiz, José Gerardo},
  journal={Computers and Mathematics with Applications},
  volume={195},
  pages={396--418},
  year={2025},
  publisher={Elsevier},
  doi={10.1016/j.camwa.2025.07.034}
}
```

<div align="center">
<i>Developed for the advancement of meshless numerical methods and scientific computing.</i><br/>
<a href="https://github.com/gstinoco/mGFD/issues">Report a Bug</a> | <a href="mailto:gerardo.tinoco@umich.mx">Contact Author</a>
</div>
