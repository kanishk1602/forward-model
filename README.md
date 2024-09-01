# MASW Surface Wave Inversion

This project focuses on implementing machine learning algorithms to enhance the accuracy of surface wave inversion, specifically for analyzing seismic data. The developed algorithms are applied to the Multi-Channel Analysis of Surface Waves (MASW) method to estimate subsurface properties such as layer velocities and depths, which are crucial for geophysical investigations in engineering and exploration.

## Project Overview

- **Algorithm Development**: Implemented machine learning algorithms to improve the accuracy of surface wave inversion. The algorithms utilize advanced mathematical techniques and data processing methods to enhance inversion modeling.
  
- **Programming**: The project is developed in C++ with a focus on high-performance computing. The Eigen library is used for dense matrix operations, optimizing the computational efficiency of the algorithms.

- **Data Processing**: Seismic surface wave data is analyzed to estimate subsurface properties. This includes calculating layer velocities and depths, which are vital for understanding subsurface structures.

- **Outcome**: The project successfully achieved high-resolution inversion results, significantly improving subsurface characterization. These results have potential applications in geophysical investigations for engineering projects and exploration activities.

## Technologies Used

- **C++**: Core programming language used for algorithm development.
- **Eigen Library**: Utilized for efficient dense matrix operations.
- **Machine Learning Algorithms**: Applied to enhance surface wave inversion accuracy.
- **Seismic Data Analysis**: Employed to estimate subsurface properties.

## Project Structure

- `src/`: Contains the source code for the project.
- `data/`: Includes sample seismic surface wave data used for testing.
- `models/`: Machine learning models developed for the inversion process.
- `results/`: Contains the output of the inversion process, including layer velocities and depth profiles.

## How to Run

1. Clone the repository:
    ```bash
    git clone https://github.com/kanishk1602/forward-model.git
    cd forward-model
    ```

2. Build the project:
    ```bash
    mkdir build
    cd build
    cmake ..
    make
    ```

3. Run the executable with your seismic data:
    ```bash
    ./MASWInversion ../data/sample_data.txt
    ```

4. The results will be saved in the `results/` directory.

## Example Output

Here is an example of the inversion results:
- High-resolution velocity profiles
- Depth estimates for subsurface layers
- Improved subsurface characterization images

## Contributing

Contributions are welcome! Feel free to open an issue or submit a pull request if you have ideas for improvements or new features.

## License

## Acknowledgments

Special thanks to the Rhythm Shah and Bharat Shekar department of Earth Sciences IIT Bombay for providing the necessary seismic data and resources for this project.

## Contact

For any queries, feel free to reach out:

- **Email**: [Kanishk Krishnan](mailto:kanishkkrishnan1602@gmail.com)
- **LinkedIn**: [kanishk1602](https://www.linkedin.com/in/kanishk1602/)
