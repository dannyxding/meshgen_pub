# MeshGen Doc

## Overview

MeshGen is a automatic command-line utility designed to generate high-quality tetrahedral meshes for ion channels. This tool integrates with third-party software like TetGen and TMSMesh2.1 to provide a seamless workflow from PQR files to finalized mesh structures suitable for computational biology studies.
## Prerequisites

- **Operating System:** Linux
- **Dependencies:** Python3
- **Third-Party Tools:** TetGen, TMSMesh2.1

Ensure that Python3 is installed on your system and that TetGen and TMSMesh2.1 are correctly installed and compiled in their respective directories.

## Installation and Setup

1. **Executable Setup:**
   Ensure the `.meshgen` executable is placed in a designated directory.

2. **PQR Files:**
   Place your original PQR files in the same directory as the `.meshgen` executable, named `<protein_id>.pqr`.

3. **Third-Party Tools:**
    - TetGen should be placed and compiled under `tool_mesh/tetgen`.
    - TMSMesh should be placed and compiled under `tool_surface/tmsmesh`.

4. **Configuration File:**
   Update the configurations in the `settings.json` file located in the same directory as the executable. This JSON file contains the settings for different protein simulations.

## Configuration Details

The `settings.json` file should contain an array of objects, each specifying the configuration for a particular protein simulation. Here is the structure of each object in the array:

- **protein_id:** Identifier for the protein (corresponds to the PQR file).
- **arguments:** Contains simulation-specific parameters:
    - **Z_1, Z_2:** Z-coordinates defining the membrane boundaries.
    - **tms_d, tms_e:** Parameters for TMSMesh.
    - **tetgen_flag:** Command-line flags for TetGen.
    - **ETA:** (if applicable) Specific parameter for the simulation.
- **extract_alg:** Algorithm used for extracting the mesh.  Possible values are "bfs", "2stepbfs" or "trail" for try both methods.

### Example Configuration Entry

```json
{
    "protein_id": "1bl8",
    "arguments": {
        "Z_1": -17,
        "Z_2": 17,
        "ETA": 20,
        "tms_d": 0.6,
        "tms_e": 0.9,
        "tetgen_flag": "-pq1.4AnQ"
    },
    "extract_alg": "bfs"
}
```

## Usage

To run the tool, use the following command:

```bash
./meshgen <protein_id>
```

Replace `<protein_id>` with the identifier of the protein you want to simulate. MeshGen will process the PQR file corresponding to the provided ID, applying the settings defined in `settings.json`.

## Output

MeshGen generates several files during execution:

- Intermediate mesh files are stored in the `off` and `poly` directories.
- The final labeled ion channel and membrane tetrahedral mesh are output to the `pvd` directory.

The execution process is cumulative, meaning MeshGen will automatically read from intermediate steps if they are consistent with the settings provided.