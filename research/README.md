# mGFD Research Laboratory :microscope:

<div align="center">

<img src="../docs/logo/logo.png" alt="mGFD logo" width="680" style="margin: 20px 0;">

**Experimental tools, benchmarks, visualizations, and reference datasets for the mGFD method.**

### :link: Quick Links
[![🧮 Model](https://img.shields.io/badge/🧮-Model-purple)](#books-mathematical-model) [![🗂️ Dataset](https://img.shields.io/badge/🗂️-Dataset-blue)](#file_cabinet-dataset-structure) [![📈 Benchmarks](https://img.shields.io/badge/📈-Benchmarks-purple)](#chart_with_upwards_trend-experiment-reproduction) [![🎬 Visualizations](https://img.shields.io/badge/🎬-Visualizations-purple)](#movie_camera-visualizations) [![👥 Team](https://img.shields.io/badge/👥-Research%20Team-blue)](#scientist-research-team) [![🏭 Partners](https://img.shields.io/badge/🏭-Industry%20Partners-0B1B3A)](#building_with_garden-collaboration-centers)

</div>

---

## :clipboard: Table of Contents
- [Mathematical Model](#books-mathematical-model)
- [Laboratory Structure](#open_file_folder-laboratory-structure)
- [Dataset Structure](#file_cabinet-dataset-structure)
- [Visualizations](#movie_camera-visualizations)
- [Experiment Reproduction](#chart_with_upwards_trend-experiment-reproduction)
- [Acknowledgements and Sponsorships](#pray-acknowledgements-and-sponsorships)
- [Research Team](#scientist-research-team)
- [Scientific References](#books-scientific-references)

---

This directory contains the entire suite of experimental tools, benchmarks, visualizations, and original datasets (real lakes) that support the scientific publication of the mGFD method.

Unlike the project root (which functions as the installable Python package), this folder (`/research/`) is the **research laboratory** where scripts consume the `mGFD` library to generate geographic studies and numerical stability analyses.

---

## :books: Mathematical Model

The mGFD method solves operator-driven PDEs on an unstructured set of nodes. A representative stationary problem is written as:

$$A u_{xx} + B u_{xy} + C u_{yy} + D u_x + E u_y + F u = -f(x,y)$$

Where boundary values are imposed through $u = \phi(x,y)$ on the boundary nodes, and the interior unknowns are calculated using local generalized finite difference stencils built from neighbor sets.

---

## :open_file_folder: Laboratory Structure

```text
research/
├── Data/                         # Point cloud datasets (by region and scale)
├── Results/                      # Generated outputs (created when running scripts)
├── batches/                      # Reproducible reference experiments
│   ├── run_Poisson.py
│   ├── run_Heat.py
│   ├── run_AdvDif.py
│   ├── run_Wave.py
│   ├── run_Perturbation.py
│   ├── run_Perturbation2.py
│   └── benchmark.py
├── benchmark_results/            # Timing and memory results
└── main.py                       # Orchestrator that runs all batches sequentially
```

---

## :file_cabinet: Dataset Structure (`Data/`)

All datasets are taken from the Author's repository [Cloud-Generation](https://github.com/gstinoco/Cloud-Generation).

### Variants (`cloud` vs `cloud_exterior`)

For each region/scale you will find two variants of point clouds:
- `*_cloud.csv`: region without holes (either because there are no islands or because the islands were filled by adding interior nodes).
- `*_cloud_exterior.csv`: region with holes (the interior islands do not have a point cloud, so they are treated as empty regions).

### Available Geographic Regions

| Region                | Abbreviation |
|-----------------------|--------------|
| `Alchichica`          | ALC          |
| `Baikal`              | BAI          |
| `Balkhash`            | BAL          |
| `Caspian`             | CAS          |
| `Catemaco`            | CAT          |
| `Chapala`             | CHA          |
| `Huron`               | HUR          |
| `Patzcuaro`           | PAT          |
| `Poopo`               | POO          |
| `Santa Maria del Oro` | SMO          |
| `Titicaca`            | TIT          |
| `Yuriria`             | YUR          |
| `Zirahuen`            | ZIR          |

---

## :movie_camera: Visualizations

### :framed_picture: Example: Lake Poopo

<div align="center">

<table>
  <tr>
    <td align="center">
      <b>Cloud</b><br/>
      <sub>Data/Poopo/&lt;scale&gt;</sub><br/><br/>
      <img src="../docs/images/Poopo_cloud.png" alt="Poopo cloud" width="260"><br/>
    </td>
    <td align="center">
      <b>Cloud exterior (with islands/holes)</b><br/>
      <sub>Data/Poopo/&lt;scale&gt; (cloud_exterior)</sub><br/><br/>
      <img src="../docs/images/Poopo_cloud_exterior.png" alt="Poopo cloud exterior" width="260"><br/>
    </td>
  </tr>
</table>

</div>

### :movie_camera: Transient Animations

<div align="center">

<table>
  <tr>
    <td align="center">
      <b>Heat (Cloud)</b><br/><br/>
      <img src="../docs/videos/Poopo_cloud_Heat.gif" alt="Heat solution (cloud, Poopo)" width="320"><br/>
      <a href="Results/Heat/Poopo/2x/cloud/"><code>Results/Heat/.../cloud/</code></a>
    </td>
    <td align="center">
      <b>Wave (Cloud exterior)</b><br/><br/>
      <img src="../docs/videos/Poopo_cloud_exterior_Wave.gif" alt="Wave solution (cloud exterior, Poopo)" width="320"><br/>
      <a href="Results/Wave/Poopo/2x/cloud_exterior/"><code>Results/Wave/.../cloud_exterior/</code></a>
    </td>
  </tr>
</table>

</div>

---

## :chart_with_upwards_trend: Experiment Reproduction

To reproduce the results generated in the scientific publications, you can run the individual batches or the main orchestrator:

```bash
# Run all batches at once
python main.py

# Run a specific experiment
python batches/run_Poisson.py
```

### Performance Benchmarks

The benchmark suite measures execution time and memory usage:

```bash
pip install pandas psutil
python batches/benchmark.py
```

The results (CSVs and text reports) are deposited in `benchmark_results/`.

---

## :scientist: Research Team

<div align="center">

### :star2: Meet the Team
*Researchers and graduate students advancing meshless computational methods*

</div>

### :busts_in_silhouette: Main Researchers

<table align="center">
  <thead>
    <tr>
      <th align="center" width="120">Photo</th>
      <th align="left">Researcher</th>
      <th align="left">Affiliation</th>
      <th align="left">Contact</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/gtinoco.webp" alt="Dr. Gerardo Tinoco Guerrero" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Dr. Gerardo Tinoco Guerrero</b> :mexico:<br/>
        <sub>Numerical Methods &amp; Computational Mathematics</sub>
      </td>
      <td>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:gerardo.tinoco@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a><br/>
        <a href="https://orcid.org/0000-0003-3119-770X"><img alt="ORCID 0000-0003-3119-770X" src="https://img.shields.io/badge/ORCID-0000--0003--3119--770X-green"></a><br/>
        <a href="https://www.researchgate.net/profile/Gerardo-Tinoco-Guerrero"><img alt="ResearchGate Profile" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </td>
    </tr>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/dmota.webp" alt="Dr. Francisco Javier Domínguez Mota" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Dr. Francisco Javier Domínguez Mota</b> :mexico:<br/>
        <sub>Applied Mathematics &amp; Finite Difference Methods</sub>
      </td>
      <td>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:francisco.mota@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a><br/>
        <a href="https://orcid.org/0000-0001-6837-172X"><img alt="ORCID 0000-0001-6837-172X" src="https://img.shields.io/badge/ORCID-0000--0001--6837--172X-green"></a><br/>
        <a href="https://www.researchgate.net/profile/Francisco-Dominguez-Mota"><img alt="ResearchGate Profile" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </td>
    </tr>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/jagt.webp" alt="Dr. José Alberto Guzmán Torres" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Dr. José Alberto Guzmán Torres</b> :mexico:<br/>
        <sub>Engineering Applications &amp; Artificial Intelligence</sub>
      </td>
      <td>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:jose.alberto.guzman@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a><br/>
        <a href="https://orcid.org/0000-0002-9309-9390"><img alt="ORCID 0000-0002-9309-9390" src="https://img.shields.io/badge/ORCID-0000--0002--9309--9390-green"></a><br/>
        <a href="https://www.researchgate.net/profile/Jose-Guzman-Torres"><img alt="ResearchGate Profile" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </td>
    </tr>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/harias.webp" alt="Dr. Heriberto Árias Rojas" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Dr. Heriberto Árias Rojas</b> :mexico:<br/>
        <sub>Engineering Applications</sub>
      </td>
      <td>
        <a href="http://www.siiia.com.mx"><img alt="Company: SIIIA MATH" src="https://img.shields.io/badge/%F0%9F%8F%A2%20Company-SIIIA%20MATH-0B1B3A"></a><br/>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:heriberto.arias@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a><br/>
        <a href="https://orcid.org/0000-0002-7641-8310"><img alt="ORCID 0000-0002-7641-8310" src="https://img.shields.io/badge/ORCID-0000--0002--7641--8310-green"></a><br/>
        <a href="https://www.researchgate.net/profile/Heriberto-Arias-Rojas"><img alt="ResearchGate Profile" src="https://img.shields.io/badge/ResearchGate-Profile-teal"></a>
      </td>
    </tr>
  </tbody>
</table>

### :mortar_board: Ph.D. Research Students

<table align="center">
  <thead>
    <tr>
      <th align="center" width="120">Photo</th>
      <th align="left">Student</th>
      <th align="left">Institution</th>
      <th align="left">Contact</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/gpj.webp" alt="Gabriela Pedraza-Jiménez" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Gabriela Pedraza-Jiménez</b><br/>
        <img alt="Ph.D. Research Student" src="https://img.shields.io/badge/Ph.D.-Research%20Student-2E8B57?style=flat-square">
      </td>
      <td>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:2220157h@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </td>
    </tr>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/eci.webp" alt="Eli Chagolla-Inzunza" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Eli Chagolla-Inzunza</b><br/>
        <img alt="Ph.D. Research Student" src="https://img.shields.io/badge/Ph.D.-Research%20Student-2E8B57?style=flat-square">
      </td>
      <td>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:1137626b@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </td>
    </tr>
  </tbody>
</table>

### :mortar_board: M.Sc. Research Students

<table align="center">
  <thead>
    <tr>
      <th align="center" width="120">Photo</th>
      <th align="left">Student</th>
      <th align="left">Institution</th>
      <th align="left">Contact</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/jlgf.webp" alt="Jorge L. González-Figueroa" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Jorge L. González-Figueroa</b><br/>
        <img alt="M.Sc. Research Student" src="https://img.shields.io/badge/M.Sc.-Research%20Student-green?style=flat-square">
      </td>
      <td>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:1718717h@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </td>
    </tr>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/cnmb.webp" alt="Christopher N. Magaña-Barocio" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Christopher N. Magaña-Barocio</b><br/>
        <img alt="M.Sc. Research Student" src="https://img.shields.io/badge/M.Sc.-Research%20Student-green?style=flat-square">
      </td>
      <td>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:1339846k@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </td>
    </tr>
  </tbody>
</table>

### :mortar_board: Undergraduate Research Students

<table align="center">
  <thead>
    <tr>
      <th align="center" width="120">Photo</th>
      <th align="left">Student</th>
      <th align="left">Institution</th>
      <th align="left">Contact</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center" width="120">
        <img src="../docs/team/mgfl.webp" alt="Maria Goretti Fraga Lopez" width="96" height="96" style="border-radius: 50%;">
      </td>
      <td>
        <b>Maria Goretti Fraga-Lopez</b><br/>
        <img alt="Undergraduate Research Student" src="https://img.shields.io/badge/Undergraduate-Research%20Student-green?style=flat-square">
      </td>
      <td>
        <a href="http://www.umich.mx"><img alt="University: UMSNH" src="https://img.shields.io/badge/%F0%9F%8E%93%20University-UMSNH-1A3A6B"></a>
      </td>
      <td>
        <a href="mailto:1702174b@umich.mx"><img alt="Contact" src="https://img.shields.io/badge/%F0%9F%93%A7-Contact-blue"></a>
      </td>
    </tr>
  </tbody>
</table>

---

## :factory: Industry Partners Supporting Innovation

<div align="center">

### :star2: Industry Partners Supporting Innovation
*Collaboration between academia and industry to accelerate real-world impact*

</div>

<div align="center">

<table align="center" width="70%">
<tr>
<td align="center">

### :factory: **SIIIA MATH**
#### *Soluciones de Ingeniería, México*

<div align="center">

[![Website](https://img.shields.io/badge/🌐-Visit%20Website-blue?style=for-the-badge)](http://www.siiia.com.mx)
[![Type](https://img.shields.io/badge/📊-R%26D%20Company-orange?style=flat-square)](http://www.siiia.com.mx)
[![Location](https://img.shields.io/badge/📍-Morelia,%20Mexico-green?style=flat-square)](http://www.siiia.com.mx)

</div>

**🎯 Focus areas:**
- Mathematical modeling & simulation
- AI/ML engineering solutions
- Technology transfer and applied R&amp;D

<div align="center">

[![Contact](https://img.shields.io/badge/📧-Partnership%20Contact-0B1B3A?style=for-the-badge)](mailto:gtinoco@siiia.com.mx)

</div>

</td>
</tr>
</table>

</div>

---

## :books: Scientific References

### :books: Core Publications (GFD / mGFD Background)

1. **Tinoco-Guerrero, G.**, Domínguez-Mota, F. J., Guzmán-Torres, J. A., & Tinoco-Ruiz, J. G. (2022). *"Numerical Solution of Diffusion Equation using a Method of Lines and Generalized Finite Differences."* **Revista Internacional de Métodos Numéricos para Cálculo y Diseño en Ingeniería**, 38(2). [DOI: 10.23967/j.rimni.2022.06.003](http://dx.doi.org/10.23967/j.rimni.2022.06.003)

### :trophy: Project Highlights

- **Meshless PDE solver core** for irregular 2D domains using Generalized Finite Differences (mGFD)
- **Reference datasets** (multiple regions, scales, and variants) enabling consistent method-to-method comparisons
- **Reproducible batch scripts** covering Poisson, Heat, Advection–Diffusion, and Wave equation families
- **Error analysis and benchmarks** with saved artifacts under `Results/` for validation and reporting

---

## :memo: Citation & License

If you use this repository in your research, please cite:

```bibtex
@software{tinoco2024mGFD,
  title={mGFD: Meshless Generalized Finite Differences for Partial Differential Equations},
  author={Tinoco-Guerrero, Gerardo and
          Domínguez-Mota, Francisco Javier and
          Guzmán-Torres, José Alberto and
          Arias-Rojas, Heriberto},
  year={2024},
  institution={Universidad Michoacana de San Nicolás de Hidalgo},
  organization={SIIIA MATH: Soluciones en ingeniería},
  url={https://github.com/gstinoco/mGFD},
  note={Meshless solver and reference datasets for irregular domains using generalized finite differences}
}
```

### :page_facing_up: License

This project is licensed under the **MIT License** (see [LICENSE](LICENSE)).

```
MIT License

Copyright (c) 2024 Gerardo Tinoco-Guerrero

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## :pray: Acknowledgments

<div align="center">

### :heart: Special Thanks
*We extend our gratitude to the institutions and partners supporting this research and open-source development*

</div>

### :classical_building: Institutional Support

<table align="center" width="100%" cellspacing="14">
  <tr>
    <td width="50%" valign="top">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 16px;">
        <div align="center">
          <b>🎓 Universidad Michoacana de San Nicolás de Hidalgo (UMSNH)</b><br/>
          <sub>Academic institution, Mexico</sub><br/><br/>
          <a href="http://www.umich.mx"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-darkred?style=flat-square"></a>
          <img alt="Type: University" src="https://img.shields.io/badge/🏷️%20Type-University-1A3A6B?style=flat-square">
          <img alt="Support: Infrastructure" src="https://img.shields.io/badge/🤝%20Support-Infrastructure-2E8B57?style=flat-square">
        </div>
        <br/>
        <b>Key support</b>
        <ul>
          <li>Academic foundation and research infrastructure</li>
          <li>Scientific training and supervision environment</li>
        </ul>
      </div>
    </td>
    <td width="50%" valign="top">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 16px;">
        <div align="center">
          <b>🏛️ Secretariat of Science, Humanities, Technology and Innovation (SECIHTI)</b><br/>
          <sub>State Secretariat, Mexico</sub><br/><br/>
          <a href="https://secihti.mx/"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-darkgreen?style=flat-square"></a>
          <img alt="Type: Government" src="https://img.shields.io/badge/🏷️%20Type-Government-2D6A4F?style=flat-square">
          <img alt="Support: Funding and Innovation" src="https://img.shields.io/badge/🤝%20Support-Funding%20%26%20Innovation-40916C?style=flat-square">
        </div>
        <br/>
        <b>Key support</b>
        <ul>
          <li>Support for science and technology initiatives</li>
          <li>Funding and innovation promotion</li>
        </ul>
      </div>
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 16px;">
        <div align="center">
          <b>🌿 Centre Internacional de Mètodes Numèrics en Enginyeria (CIMNE)</b><br/>
          <sub>Industry, Spain</sub><br/><br/>
          <a href="https://aulas.cimne.com/aula/aula-morelia/"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-orange?style=flat-square"></a>
          <img alt="Type: Research Center" src="https://img.shields.io/badge/🏷️%20Type-Research%20Center-EE9B00?style=flat-square">
          <img alt="Support: Collaboration" src="https://img.shields.io/badge/🤝%20Support-Collaboration-CA6702?style=flat-square">
        </div>
        <br/>
        <b>Key support</b>
        <ul>
          <li>International collaboration in numerical methods</li>
          <li>Computational engineering research environment</li>
        </ul>
      </div>
    </td>
    <td width="50%" valign="top">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 16px;">
        <div align="center">
          <b>🏭 SIIIA MATH: Soluciones en Ingeniería</b><br/>
          <sub>Industry, México</sub><br/><br/>
          <a href="http://www.siiia.com.mx"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-blue?style=flat-square"></a>
          <img alt="Type: Industry Partner" src="https://img.shields.io/badge/🏷️%20Type-Industry%20Partner-0B1B3A?style=flat-square">
          <img alt="Support: Technology Transfer" src="https://img.shields.io/badge/🤝%20Support-Technology%20Transfer-1D3557?style=flat-square">
        </div>
        <br/>
        <b>Key support</b>
        <ul>
          <li>Industry-driven applied research and development</li>
          <li>Technology transfer and practical engineering impact</li>
        </ul>
      </div>
    </td>
  </tr>
</table>

### :building_with_garden: Research Centers & Collaborations

<div align="center">

<table align="center" width="100%" cellspacing="14">
  <tr>
    <td width="50%" valign="top">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 16px;">
        <div align="center">
          <b>🌿 Aula CIMNE-Morelia</b><br/>
          <sub>Research collaboration space</sub><br/><br/>
          <a href="https://aulas.cimne.com/aula/aula-morelia/"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-orange?style=flat-square"></a>
          <img alt="Area: Numerical Methods" src="https://img.shields.io/badge/🧮%20Area-Numerical%20Methods-EE9B00?style=flat-square">
          <img alt="Collaboration: Applied Computing" src="https://img.shields.io/badge/🤝%20Collaboration-Applied%20Computing-CA6702?style=flat-square">
        </div>
        <br/>
        <b>Collaboration highlights</b>
        <ul>
          <li>Numerical methods and computational engineering environment</li>
          <li>Academic–industry collaboration and training activities</li>
        </ul>
      </div>
    </td>
    <td width="50%" valign="top">
      <div style="border: 1px solid #d0d7de; border-radius: 12px; padding: 16px;">
        <div align="center">
          <b>🎓 UMSNH</b><br/>
          <sub>Academic collaboration</sub><br/><br/>
          <a href="http://www.umich.mx"><img alt="Website" src="https://img.shields.io/badge/🌐-Website-darkred?style=flat-square"></a>
          <img alt="Type: University" src="https://img.shields.io/badge/🏷️%20Type-University-1A3A6B?style=flat-square">
          <img alt="Support: Research Infrastructure" src="https://img.shields.io/badge/🤝%20Support-Research%20Infrastructure-2E8B57?style=flat-square">
        </div>
        <br/>
        <b>Collaboration highlights</b>
        <ul>
          <li>Institutional infrastructure supporting research and training</li>
          <li>Graduate formation and supervision for scientific computing</li>
        </ul>
      </div>
    </td>
  </tr>
</table>

</div>

### :computer: Technology Communities

<div align="center">

| :package: Framework | :busts_in_silhouette: Community | :star: Contribution |
|:---:|:---:|:---:|
| [![NumPy](https://img.shields.io/badge/NumPy-Scientific%20Computing-013243?style=flat-square&logo=numpy)](https://numpy.org/) | **NumPy Community** | Array computing foundation |
| [![SciPy](https://img.shields.io/badge/SciPy-Scientific%20Computing-8CAAE6?style=flat-square&logo=scipy)](https://scipy.org/) | **SciPy Community** | Numerical algorithms |
| [![Matplotlib](https://img.shields.io/badge/Matplotlib-Visualization-11557C?style=flat-square)](https://matplotlib.org/) | **Matplotlib Community** | Scientific visualization |
| [![Joblib](https://img.shields.io/badge/Joblib-Parallelism-2E8B57?style=flat-square)](https://joblib.readthedocs.io/) | **Joblib Community** | Parallel computing utilities |

</div>

---



## :books: Scientific References

1. **Tinoco-Guerrero, G.**, Domínguez-Mota, F. J., Guzmán-Torres, J. A., & Tinoco-Ruiz, J. G. (2022). *"Numerical Solution of Diffusion Equation using a Method of Lines and Generalized Finite Differences."* **Revista Internacional de Métodos Numéricos para Cálculo y Diseño en Ingeniería**, 38(2). [DOI: 10.23967/j.rimni.2022.06.003](http://dx.doi.org/10.23967/j.rimni.2022.06.003)

If you use this laboratory in your research, please cite:

```bibtex
@software{tinoco2024mGFD,
  title={mGFD: Meshless Generalized Finite Differences for Partial Differential Equations},
  author={Tinoco-Guerrero, Gerardo and Domínguez-Mota, Francisco Javier and Guzmán-Torres, José Alberto and Arias-Rojas, Heriberto},
  year={2024},
  institution={Universidad Michoacana de San Nicolás de Hidalgo},
  organization={SIIIA MATH: Soluciones en ingeniería},
  url={https://github.com/gstinoco/mGFD},
  note={Meshless solver and reference datasets for irregular domains using generalized finite differences}
}
```

<div align="center">
<i>If this scientific laboratory helps in your research, consider giving a star ⭐️ to the main project!</i><br/><br/>
<i>Built with ❤️ by the scientific computing community in Morelia, México.</i><br/>
<a href="https://github.com/gstinoco/mGFD/issues">Report a Bug</a> | <a href="mailto:gerardo.tinoco@umich.mx">Contact Author</a>
</div>
