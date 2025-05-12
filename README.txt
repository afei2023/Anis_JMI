Anisotropic(VTI) Joint Migration Inversion (JMI)
﻿
This project demonstrates a basic workflow of the Joint Migration Inversion (JMI) algorithm, consisting of three main steps. Please follow the order of execution (a → b → c) to reproduce the example results.
﻿
Step a: Build the Salt-Dome Model
- Construct a synthetic salt-dome velocity model.
- This model serves as the initial true model for subsequent simulations and inversions.
﻿
Step b: Perform Wavefield Simulation
- Simulate seismic wavefields using the constructed model.
- Both primary reflections (first-order waves) and interbed multiples (higher-order multiples) are generated to provide more complete illumination of the subsurface.
﻿
Step c: Conduct Joint Migration Inversion
- Perform reflectivity migration and velocity model inversion in a closed-loop manner.
- Update the reflectivity and velocity iteratively to improve imaging accuracy and model fidelity.
﻿
Notes:
- Please ensure all dependencies are installed before running the scripts.
- Execution order is critical to maintain data consistency across steps.
- For detailed parameter settings and further customization, refer to comments within each script.
﻿
We hope this example helps you understand and apply the JMI workflow efficiently.