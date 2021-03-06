# Li-ion-Battery-Simulation-with-MATLAB
This is a numerical solver based on finite difference method for pseudo two-dimensional model. 

And the folder [/images](./images) contains a demo (the input parameters are *not real cases*, you may treat it as a **toy example**).

Note that it currently has following **limitations** to be improved:
- work for only *small rates* due to used simplifications
- consider only *isothermal* conditions 
- has *not verify* exactness

References (more strict/convincing methods/results):

1. [Modeling of Galvanostatic Charge and Discharge of the Lithium/Polymer/Insertion Cell.](http://jes.ecsdl.org/content/140/6/1526.abstract)
1. [State-of-Charge Estimation for Lithium-Ion Batteries Using Neural Networks and EKF.](https://ieeexplore.ieee.org/abstract/document/5416277)
1. [Lionsimba: a Matlab Framework Based on a Finite Volume Model Suitable for Li-ion Battery Design, Simulation, and Control.](https://github.com/lionsimbatoolbox/LIONSIMBA)
1. [基于电化学模型的锂离子电池多尺度建模及其简化方法.](http://wulixb.iphy.ac.cn/CN/abstract/abstract71145.shtml)
1. [Lithium-ion Battery Thermal-electrochemical Model-based State Estimation Using Orthogonal Collocation and a Modified Extended Kalman Filter.](https://arxiv.org/pdf/1506.08689.pdf)
1. [Coordinate Transformation, Orthogonal Collocation, Model Reformulation and Simulation of Electrochemical-Thermal Behavior of Lithium-Ion Battery Stacks.](http://jes.ecsdl.org/content/158/12/A1461.full.pdf?casa_token=x4z_ND9r7HUAAAAA:3kQ8V2vdedQcDfRfS3qa1OqE6k5JU2AtRwomVQa3rEP48JLtPzbmVmXqDnQLk3TyjllH_w-fbZVMMg)
