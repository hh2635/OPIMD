### 1. Configuration
Navigate to the parameter file to configure your simulation settings:
```bash
# Locate the parameter file
# Note: Edit this file with any text editor (vim, nano, or vs code)
vi src/inputs/parameters.json

### 2. Simulation
cd src/auto/
chmod +x submit_job.sh  # Ensure the script is executable
./submit_job.sh
