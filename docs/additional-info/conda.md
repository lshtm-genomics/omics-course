
# Conda Environments: Why and How to Use Them

## What is a Conda Environment?

Conda is an open-source package manager and environment manager that allows users to install multiple versions of software packages and their dependencies in isolated environments. A conda environment is essentially a directory that contains a specific collection of conda packages.

## Why Use Conda Environments?

1. **Isolation**: Multiple projects may require different versions of libraries or even the Python interpreter itself. By using separate environments, you can ensure that each project has its own dependencies, isolated from each other.
2. **Avoid Conflicts**: Dependencies of one project might conflict with another. Environments ensure that such conflicts are avoided by isolating dependencies.
3. **Reproducibility**: Sharing a project with someone else (or even with your future self) becomes much simpler when you can provide a list of exact package versions used in your environment.
4. **Flexibility**: You can easily switch between different versions of a package or Python itself, depending on the project's needs.
5. **Clean and Safe Experimentation**: If you want to test out a new package or update an existing one, you can do so in a new environment. If something goes wrong, your main working environment remains unaffected.

## How Does It Work?

Conda environments work by creating isolated spaces, each with its own installation directories. This allows you to have different versions of a package (or even Python itself) in different environments. When you activate an environment, conda adjusts your `PATH` so that the selected environment's executables are used.

## Some Basic Conda Commands:

1. **Creating a New Environment**:  
   `conda create --name myenv python=3.7`
2. **Activating an Environment**:  
   - On macOS and Linux: `source activate myenv` or `conda activate myenv`
3. **Deactivating an Environment**:  
   `conda deactivate`
4. **Listing All Environments**:  
   `conda env list` or `conda info --envs`
5. **Installing a Package in an Active Environment**:  
   `conda install numpy`
6. **Exporting an Environment to a YAML File**:  
   `conda env export > environment.yml`
7. **Creating an Environment from a YAML File**:  
   `conda env create -f environment.yml`
8. **Removing an Environment**:  
   `conda env remove --name myenv`

## Environments for the course

Click [here](https://github.com/jodyphelan/genomics_course/tree/master/conda_env) for links to environment file used for the practicals.

## Conclusion

Conda environments are essential tools for data scientists and developers to manage dependencies and ensure project reproducibility. By understanding and utilizing them, you can maintain a more organized, conflict-free, and efficient development workflow.
