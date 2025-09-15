# mGFD CloudGenerator 2.0 :cloud:

<div align="center">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-black.svg)](https://github.com/gstinoco/CloudGen) [![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/) [![Flask](https://img.shields.io/badge/Flask-2.3+-green.svg)](https://flask.palletsprojects.com/) [![OpenCV](https://img.shields.io/badge/OpenCV-4.8+-red.svg)](https://opencv.org/) [![GMSH](https://img.shields.io/badge/GMSH-4.11+-orange.svg)](https://gmsh.info/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Advanced Web Platform for Generating Unstructured Clouds of Points**

*Complete solution for meshless Generalized Finite Difference Method (mGFD) applications*

### :link: Quick Links
[![:rocket: Quick Start](https://img.shields.io/badge/🚀-Quick%20Start-green)](#rocket-quick-start) [![:computer: Features](https://img.shields.io/badge/💻-Features-blue)](#sparkles-features) [![:busts_in_silhouette: Team](https://img.shields.io/badge/👥-Research%20Team-blue)](#scientist-research-team)

</div>

---

## :clipboard: Table of Contents
- [Overview](#star2-overview)
- [Features](#sparkles-features)
- [Installation & Setup](#package-installation--setup)
- [Quick Start](#rocket-quick-start)
- [Usage Guide](#book-usage-guide)
- [API Documentation](#gear-api-documentation)
- [Examples](#bulb-examples)
- [Project Architecture](#open_file_folder-project-architecture)
- [Configuration](#wrench-configuration)
- [Scientific Background](#books-scientific-background)
- [Contributing](#handshake-contributing)
- [Research Team](#scientist-research-team)
- [Citation & License](#memo-citation--license)
- [Contact](#email-contact--support)

---

## :star2: Overview

The **mGFD CloudGenerator 2.0** is a comprehensive web-based platform designed for generating optimized unstructured clouds of points specifically tailored for the meshless Generalized Finite Difference Method (mGFD). This advanced tool combines interactive image processing capabilities with sophisticated cloud generation algorithms to provide researchers and engineers with a complete solution for numerical simulations.

### :gear: Key Capabilities
- **:art: Interactive Contour Creation**: Advanced image segmentation with click-based region detection
- **:cloud: Optimized Cloud Generation**: High-quality point cloud generation using GMSH integration
- **:chart_with_upwards_trend: Real-time Visualization**: Interactive canvas with zoom, pan, and multi-region support
- **:floppy_disk: Multiple Export Formats**: CSV data export with PNG/SVG visualizations
- **:globe_with_meridians: Web-based Interface**: Modern, responsive design accessible from any browser

### :microscope: Applications

| Field | Application | Use Case |
|-------|-------------|----------|
| **Computational Fluid Dynamics** :ocean: | Flow Simulation | Irregular domain discretization, boundary layer modeling |
| **Structural Engineering** :building_construction: | Stress Analysis | Complex geometry meshing, crack propagation studies |
| **Heat Transfer** :fire: | Thermal Analysis | Non-uniform domain discretization, interface problems |
| **Environmental Modeling** :herb: | Pollution Transport | Irregular terrain modeling, contaminant dispersion |
| **Biomedical Engineering** :microscope: | Tissue Modeling | Organ geometry discretization, drug delivery simulation |

---

## :sparkles: Features

### :art: ContourCreator Module
- **Interactive Image Processing**: Upload and process images (PNG, JPG, JPEG, GIF, BMP)
- **Advanced Segmentation**: OpenCV-based region detection with customizable tolerance
- **Multi-region Support**: Detect and manage multiple regions with color-coded visualization
- **Canvas Manipulation**: Zoom, pan, and precise click-based region selection
- **Real-time Preview**: Instant visualization of detected contours and regions

### :cloud: CloudGenerator Module
- **GMSH Integration**: Professional mesh generation with adaptive sizing algorithms
- **Node Classification**: Automatic classification (interior, boundary, interface nodes)
- **Memory Optimization**: Efficient processing for large datasets with garbage collection
- **Asynchronous Processing**: Background cloud generation with real-time status updates
- **Quality Visualization**: High-resolution PNG and scalable SVG output formats

### :globe_with_meridians: Web Interface
- **Modern Design**: Responsive interface with glassmorphism effects and smooth animations
- **Drag & Drop**: Intuitive file upload with progress indicators
- **Real-time Feedback**: Live status updates and error handling
- **Professional Logging**: Comprehensive logging system with file rotation
- **Cross-platform**: Compatible with all modern web browsers

---

## :package: Installation & Setup

### :computer: System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **Python** | 3.8+ | 3.9+ |
| **RAM** | 4 GB | 8 GB+ |
| **CPU** | 2 cores | 4+ cores |
| **Storage** | 1 GB | 5 GB+ (for datasets) |
| **OS** | Windows/Linux/macOS | Linux (optimal performance) |

### :package: Dependencies

The project uses the following main dependencies:

```python
# Core web framework
Flask >= 2.3.0              # Web application framework
Werkzeug >= 2.3.0           # WSGI utilities

# Computer vision and image processing
opencv-python >= 4.8.0      # Image processing and segmentation
Pillow >= 10.0.0             # Image manipulation

# Scientific computing
numpy >= 1.24.0             # Numerical computations
pandas >= 2.0.0              # Data manipulation
scipy >= 1.10.0              # Scientific algorithms

# Mesh generation
gmsh >= 4.11.0               # Geometric modeling and mesh generation

# Visualization
matplotlib >= 3.7.0         # Scientific plotting
seaborn >= 0.12.0            # Statistical visualization
```

### :wrench: Installation Steps

#### Method 1: Direct Installation
```bash
# Clone the repository
git clone https://github.com/gstinoco/CloudGen.git
cd CloudGen

# Install dependencies
pip install -r requirements.txt
```

#### Method 2: Virtual Environment (Recommended)
```bash
# Create virtual environment
python -m venv mGFD_env
source mGFD_env/bin/activate  # On Windows: mGFD_env\Scripts\activate

# Clone and install
git clone https://github.com/gstinoco/CloudGen.git
cd CloudGen
pip install -r requirements.txt
```

#### Method 3: Conda Environment
```bash
# Create conda environment
conda create -n mGFD_cloudgen python=3.9
conda activate mGFD_cloudgen

# Clone and install
git clone https://github.com/gstinoco/CloudGen.git
cd CloudGen
pip install -r requirements.txt
```

### :white_check_mark: Installation Verification

```bash
# Test installation
python -c "import flask, cv2, gmsh, numpy, pandas; print('✅ Installation successful!')"

# Run the application
python app.py
```

The application will be available at `http://localhost:8080`

---

## :rocket: Quick Start

### 1. :arrow_forward: Launch the Application

```bash
# Navigate to project directory
cd mGFD-CloudGenerator

# Start the web server
python app.py
```

### 2. :art: Create Contours (ContourCreator)

1. **Upload Image**: Navigate to ContourCreator and upload your image
2. **Detect Regions**: Click on regions of interest to detect contours
3. **Manage Regions**: Add, remove, or modify detected regions
4. **Export Data**: Save coordinates as CSV files

### 3. :cloud: Generate Point Clouds (CloudGenerator)

1. **Upload CSV**: Use the CSV file from ContourCreator or upload your own
2. **Configure Parameters**: Set generation options (regions inside/outside)
3. **Generate Cloud**: Start the cloud generation process
4. **Download Results**: Get CSV data and visualization files

### 4. :chart_with_upwards_trend: Analyze Results

- **CSV Files**: Node coordinates with region classification
- **PNG Images**: High-resolution visualizations
- **SVG Files**: Scalable vector graphics for publications

---

## :book: Usage Guide

### ContourCreator Workflow

#### Step 1: Image Upload
```javascript
// Supported formats: PNG, JPG
// Maximum file size: 10MB
// Drag & drop or click to browse
```

#### Step 2: Region Detection
- Click on image regions to detect contours
- Adjust tolerance for segmentation sensitivity
- Use zoom and pan for precise selection
- Multiple regions supported with color coding

#### Step 3: Region Management
- **Add Region**: Click "Add Region" after detection
- **Toggle Visibility**: Show/hide regions individually
- **Delete Region**: Remove unwanted regions
- **Clear All**: Reset all detected regions

#### Step 4: Data Export
- **Single Region**: Export individual region coordinates
- **All Regions**: Export complete dataset with region labels
- **CSV Format**: Normalized coordinates (0-1 range)

### CloudGenerator Workflow

#### Step 1: CSV Upload
```csv
# Expected CSV format:
x,y,region
0.1,0.2,1
0.3,0.4,1
0.5,0.6,2
```

#### Step 2: Configuration
- **Regions Inside**: Generate points inside detected regions
- **Regions Outside**: Generate points outside detected regions
- **Adaptive Sizing**: Automatic point density optimization

#### Step 3: Generation Process
- **Asynchronous Processing**: Background generation with status updates
- **Memory Management**: Optimized for large datasets
- **Error Handling**: Comprehensive error reporting

#### Step 4: Results
- **Node Classification**: Interior, boundary, and interface nodes
- **Multiple Formats**: CSV data with PNG/SVG visualizations
- **Quality Metrics**: Point distribution analysis

---

## :open_file_folder: Project Architecture

```
mGFD-CloudGenerator/
├── 📄 app.py                    # Main Flask application
├── 📄 cloud_generation.py       # Cloud generation algorithms
├── 📄 reduce_points.py          # Point reduction utilities
├── 📄 requirements.txt          # Python dependencies
├── 📁 templates/               # HTML templates
│   ├── 🏠 home.html            # Landing page
│   ├── 🎨 contour_creator.html # ContourCreator interface
│   ├── ☁️ cloud_generator.html  # CloudGenerator interface
│   └── ℹ️ about.html            # About page
├── 📁 static/                  # Static assets
│   ├── 🎨 css/styles.css       # Main stylesheet (6000+ lines)
│   ├── 📜 js/                  # JavaScript modules
│   │   ├── contour_creator.js  # ContourCreator functionality
│   │   ├── cloud_generator.js  # CloudGenerator functionality
│   │   └── navbar.js           # Navigation components
│   ├── 🖼️ images/              # Logos and assets
│   └── 📊 examples/            # Sample data files
├── 📁 uploads/                 # Temporary file storage
├── 📁 output/                  # Generated results
└── 📁 logs/                    # Application logs
```

### Core Modules

#### Flask Application (`app.py`)
- **Web Framework**: Flask-based REST API
- **File Management**: Upload handling and cleanup
- **Image Processing**: OpenCV integration for segmentation
- **Async Processing**: Background task management
- **Logging System**: Professional logging with rotation

#### Cloud Generation (`cloud_generation.py`)
- **GMSH Integration**: Geometric modeling and mesh generation
- **Adaptive Algorithms**: Point density optimization
- **Node Classification**: Interior/boundary/interface detection
- **Memory Management**: Efficient processing for large datasets
- **Visualization**: High-quality PNG/SVG output

#### Frontend (`static/`)
- **Modern UI**: Responsive design with CSS Grid/Flexbox
- **Interactive Canvas**: HTML5 Canvas with zoom/pan capabilities
- **Real-time Updates**: WebSocket-like status monitoring
- **File Handling**: Drag & drop with progress indicators

---

## :handshake: Contributing

We welcome contributions from the research community! Here's how you can help:

### :bug: Bug Reports
1. **Search Existing Issues**: Check if the bug has been reported
2. **Create Detailed Report**: Include steps to reproduce, expected vs actual behavior
3. **Provide Context**: Operating system, Python version, browser details
4. **Include Logs**: Attach relevant log files from the `logs/` directory

### :bulb: Feature Requests
1. **Describe the Feature**: Clear description of the proposed functionality
2. **Justify the Need**: Explain how it benefits the research community
3. **Provide Examples**: Include use cases and expected behavior
4. **Consider Implementation**: Suggest possible approaches if applicable

### :computer: Code Contributions

#### Development Setup
```bash
# Fork the repository
git clone https://github.com/yourusername/CloudGen.git
cd CloudGen

# Create development environment
python -m venv dev_env
source dev_env/bin/activate
pip install -r requirements.txt

# Create feature branch
git checkout -b feature/your-feature-name
```

#### Coding Standards
- **Python Style**: Follow PEP 8 guidelines
- **Documentation**: Include comprehensive docstrings
- **Testing**: Add unit tests for new functionality
- **Logging**: Use the existing logging framework
- **Error Handling**: Implement robust error handling

#### Pull Request Process
1. **Update Documentation**: Ensure README and docstrings are current
2. **Test Thoroughly**: Verify functionality across different scenarios
3. **Follow Conventions**: Maintain consistent code style
4. **Describe Changes**: Provide clear PR description with examples

### :memo: Documentation
- **API Documentation**: Help improve endpoint documentation
- **User Guides**: Create tutorials and usage examples
- **Scientific Papers**: Contribute to research publications
- **Translations**: Help translate documentation to other languages

---

## :busts_in_silhouette: Research Team

### :microscope: **Principal Researchers**

#### :man_office_worker: **Dr. Gerardo Tinoco-Guerrero** - *Principal Researcher*
> :mortar_board: **Ph.D. in Physical Engineering Sciences** | Universidad Michoacana de San Nicolás de Hidalgo

- :bar_chart: **Leadership**: Project coordination and scientific direction
- :email: **Contact**: [gerardo.tinoco@umich.mx](mailto:gerardo.tinoco@umich.mx)

#### :man_scientist: **Dr. José Alberto Guzmán-Torres** - *Co-Researcher*
> :mortar_board: **Ph.D. in Physical Engineering Sciences** | Universidad Michoacana de San Nicolás de Hidalgo

- :bulb: **Contribution**: Technical implementation and validation
- :email: **Contact**: [jose.alberto.guzman@umich.mx](mailto:jose.alberto.guzman@umich.mx)

#### :man_scientist: **Dr. Francisco Javier Domínguez-Mota** - *Co-Researcher*
> :mortar_board: **Ph.D. in Mathematical Sciences** | Universidad Michoacana de San Nicolás de Hidalgo

- :abacus: **Contribution**: Mathematical rigor and theoretical analysis
- :email: **Contact**: [francisco.mota@umich.mx](mailto:francisco.mota@umich.mx)

---

### :mortar_board: **Graduate Students**

<div align="center">

| :woman_student::man_student: **Student** | :chart_with_upwards_trend: **Status** |
|:---|:---:|
| **Gabriela Pedraza-Jiménez** | ![PhD](https://img.shields.io/badge/Ph.D.-Candidate-purple) |
| **Eli Chagolla-Inzunza** | ![PhD](https://img.shields.io/badge/Ph.D.-Candidate-purple) |
| **Ángel E. Calvillo-Vázquez** | ![MSc](https://img.shields.io/badge/M.Sc.-Student-green) |
| **Jorge L. González-Figueroa** | ![MSc](https://img.shields.io/badge/M.Sc.-Student-green) |
| **Christopher N. Magaña-Barocio** | ![MSc](https://img.shields.io/badge/M.Sc.-Student-green) |

</div>

### :star2: **Team Contributions**

- :microscope: **Ph.D. Candidates**: Advanced research in computational methods and numerical analysis
- :man_technologist: **M.Sc. Students**: Software development, algorithm implementation, and data processing
- :handshake: **Collaborative Approach**: Interdisciplinary team combining mathematics, engineering, and computer science
- :books: **Academic Excellence**: All members affiliated with Universidad Michoacana de San Nicolás de Hidalgo

### :handshake: Sponsors

<div align="center">

| :office: **Organization** | :globe_with_meridians: **Type** | :round_pushpin: **Location** | :link: **Website** |
|:---:|:---:|:---:|:---:|
| **SIIIA MATH** | R&D Company | Morelia, México | [![Website](https://img.shields.io/badge/🌐-Visit%20Site-blue)](http://siiia.com.mx/) |
| **UMSNH** | Public University | Morelia, México | [![Website](https://img.shields.io/badge/🌐-Visit%20Site-green)](https://umich.mx/) |
| **SECIHTI** | Government Agency | México | [![Website](https://img.shields.io/badge/🌐-Visit%20Site-orange)](https://secihti.mx/) |

</div>

---

#### :factory: **[SIIIA MATH: Soluciones de Ingeniería](https://siiia.com.mx/)**
> *Transforming complex challenges into innovative solutions through advanced mathematical models and cutting-edge technology*

- :dart: **Specialization**: Artificial Intelligence & Engineering Solutions
- :trophy: **Experience**: 12+ years in mathematical modeling
- :rocket: **Projects**: 15+ completed projects with satisfied clients

#### :mortar_board: **[Universidad Michoacana de San Nicolás de Hidalgo (UMSNH)](https://umich.mx/)**
> *Leading public university fostering scientific research and academic excellence*

- :books: **Founded**: 1917 - Over 100 years of academic tradition
- :microscope: **Research**: Advanced computational mathematics and engineering
- :busts_in_silhouette: **Community**: Home to our research team and scientific development

#### :classical_building: **[Secretaría de Ciencia, Humanidades, Tecnología e Innovación (SECIHTI)](https://secihti.mx/)**
> *Mexican government agency promoting science, technology and innovation for national development*

- :mexico: **Mission**: Advancing Mexico's scientific and technological capabilities
- :bulb: **Focus**: Supporting research projects and innovation initiatives
- :star2: **Impact**: Fostering collaboration between academia and industry

### Publications
*No related publications yet. Research is ongoing and publications are in preparation.*

---

## :memo: Citation & License

### Citation

If you use mGFD CloudGenerator in your research, please cite:

```bibtex
@software{tinoco2025mGFD,
  title={mGFD CloudGenerator 2.0: Advanced Web Platform for Generating Unstructured Clouds of Points},
  author={Tinoco-Guerrero, Gerardo},
  year={2025},
  url={https://github.com/gstinoco/CloudGen},
  version={2.0}
}
```

### License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2025 Gerardo Tinoco-Guerrero

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

### Acknowledgments

- **GMSH Development Team**: For the excellent mesh generation library
- **OpenCV Community**: For computer vision and image processing tools
- **Flask Development Team**: For the lightweight web framework
- **Scientific Python Community**: For NumPy, SciPy, and Matplotlib
- **Research Community**: For feedback and contributions

---

## :email: Contact & Support

### :mortar_board: Academic Inquiries

**Research Collaboration**
- **Email**: gerardo.tinoco@umich.mx
- **Institution**: Universidad Michoacana de San Nicolás de Hidalgo
- **Topics**: Meshless methods, numerical analysis, scientific computing

### :computer: Technical Support

**Bug Reports & Feature Requests**
- **GitHub Issues**: [Create an issue](https://github.com/gstinoco/CloudGen/issues)
- **Documentation**: Check the [Wiki](https://github.com/gstinoco/CloudGen/wiki)
- **Discussions**: Join [GitHub Discussions](https://github.com/gstinoco/CloudGen/discussions)

### :handshake: Community

**Stay Connected**
- **GitHub**: Follow the repository for updates
- **Research Gate**: Connect with the research team
- **Academic Networks**: Find us on academic social platforms

### :question: FAQ

**Common Questions**

**Q: What file formats are supported for images?**
A: PNG, JPG, JPEG, GIF, and BMP formats up to 10MB.

**Q: Can I use this for commercial applications?**
A: Yes, the MIT license allows commercial use with proper attribution.

**Q: How do I cite this work in my research?**
A: Use the BibTeX citation provided in the Citation section.

**Q: Is there a limit on the number of points generated?**
A: The limit depends on your system memory. The tool is optimized for large datasets.

**Q: Can I contribute new algorithms?**
A: Absolutely! We welcome contributions. Please see the Contributing section.

---

<div align="center">

**Made with :heart: for the Scientific Computing Community**

*Advancing meshless methods through open-source collaboration*

[![GitHub stars](https://img.shields.io/github/stars/gstinoco/CloudGen?style=social)](https://github.com/gstinoco/CloudGen/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/gstinoco/CloudGen?style=social)](https://github.com/gstinoco/CloudGen/network/members)
[![GitHub watchers](https://img.shields.io/github/watchers/gstinoco/CloudGen?style=social)](https://github.com/gstinoco/CloudGen/watchers)

</div>