# MolVR

MolVR is a Python-based application designed to visualize molecular structures in 3D, offering an immersive experience by supporting Head-Mounted Display (HMD) technology. It allows users to load molecular data from `.mol` files and view them in a virtual reality environment, making it a powerful tool for education, research, and development in fields like chemistry, biochemistry, and molecular biology.

## Features

- **3D Visualization**: Render molecular structures in 3D from `.mol` files.
- **VR Support**: Experience molecular models in virtual reality for an immersive learning and research tool.
- **Cross-Platform**: Built with Python, making it relatively easy to run on various operating systems.
- **Extensible**: Open source and designed for easy extension and customization.

## Installation

Before installing MolVR, ensure you have Python 3.8 or later installed on your system. MolVR also requires specific Python libraries and a compatible VR SDK.

1. **Clone the Repository**

   ```bash
   git clone git@github.com:lukebellis/Molvr.git
   cd Molvr
Set Up a Virtual Environment (Optional but recommended)


python -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
Install Dependencies


pip install -r requirements.txt

## Usage

To start visualizing molecular structures, you can run the main script from the project directory:

python src/main.py
Replace main.py with the appropriate script name if your entry point differs.

## Contributing
Contributions to MolVR are welcome! If you're interested in helping to improve the project, please feel free to fork the repository, make your changes, and submit a pull request.

## License
MolVR is open source and distributed under the MIT license. See the LICENSE file for more details.

## Acknowledgments
Thanks to the RDKit and PyOpenGL communities for providing the essential libraries used in this project.
Special thanks to all contributors and supporters of the MolVR project.

