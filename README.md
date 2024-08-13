# Computational-Model-for-Li-SVO-CFx-Batteries-used-in-Implantable-Medical-Devices
This regulatory science tool (RST) is a computational model for predicting implantable lithium battery temperature, remaining capacity and longevity.

# Technical Description:
This tool is a coupled electro-thermal dynamic model that estimates performance of a lithium sulfur vanadium oxide - carbon monofluoride (Li/SVO-CFx) hybrid cathode primary battery over time. This lithium battery chemistry is particularly prevalent in implantable cardiac devices such as implantable cardioverter defibrillators (ICDs). This tool simplifies assessing how device electrical load impacts the battery lifespan. The input to this model is discharge current profile as a function of time and the output is battery terminal voltage. The battery characteristic Open Circuit Voltage-State of charge (OCV-SOC) curve is also provided to the model. This tool predicts three internal states: battery transient voltage, depth of discharge (DOD), and internal temperature. Additionally, the battery terminal voltage is the model output. The user can adjust parameters such as battery nominal capacity, initial state of charge, simulation duration, and the load current profile. Acceptable parameters and input ranges are provided in both the “User Manual” document and the code.

# Intended Purpose:
This model is intended to be used for predicting implantable battery performance parameters including the remaining capacity, terminal voltage, and surface temperature of implantable battery during its operation under various temperature and electrical load conditions. For example, in an Implantable Cardioverter Defibrillator (ICD), the model simulates battery performance and safety under specific operating conditions for ICDs, including the impact of high-energy shocks (defibrillation) and low-energy device power supply (housekeeping). The simulation tool enables enhanced battery safety and accurate end-of-life prediction, thus extending the safety and reliability of implantable medical devices and reducing the need for invasive replacement surgery. The tool output might be used to informed adjustment of the battery setting prior to device implementation and also the real-time device operation.

# Testing:
We investigate the performance of our model by comparing the predicted voltage-capacity curve with experimental bench test data from the literature. The Mean Absolute Percentage Error (MAPE) between the predicted and observed voltage values across the two experimental conditions (30 µA and 188 µA) is approximately 6.35%. This error can be significantly reduced by optimizing the model parameters for the specific battery used in the experiment. As a use case we simulated the Li/SVO-CFx battery dynamics for an actual ICD load parameters taken from the literature. Moreover, we quantified the model uncertainty on the battery initial state estimation by conducting a Monte Carlo simulation and Cram'er-Rao lower bound. The linear least-squares initial states estimator could estimate the initial battery DOD and initial transient voltage with an error of ±5% and ±0.02% respectively.

# Limitations:
1. The model has only been specifically tested and evaluated for Li/SVO-CFx battery chemistry, not other chemistries.
2. The model has only been evaluated using simulated load conditions representing typical ICD usage, not other medical devices.
3. The model assumes uniform battery temperature, constant material properties, and that the battery’s chemical reaction is much faster than its electrical response.

# Supporting Documentation:
The link to the RST’s journal publication:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11218594/
