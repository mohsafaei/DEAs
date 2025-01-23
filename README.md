# DEAs
## [In silico actuation performance investigation of dielectric elastomers with TPMS geometries](https://doi.org/10.1016/j.euromechsol.2024.105540)


Part of the code for creating TPMS structures was developed by [Fayyaz Nosouhi](dehnavifn@gmail.com) and [Saeed Khaleghi](saeedkhaleghi123@gmail.com).  
The remaining parts of the code, including the definition of DEA composites and their implementation in ABAQUS, were developed by Mohammad Ali Safaei (mohammadsf1998@gmail.com).<br />
This Python script requires the use of a UEL subroutine, which implements the constitutive equations for a DEA element. In this regard, the UEL subroutine developed by Ehsan Hajiesmaili was utilized. 
The UEL file and its explanation are available in the following research article: [Link](https://pubs.aip.org/aip/jap/article/129/15/151102/1025587/Dielectric-elastomer-actuators)
University of Tehran, December 2024





## Steps to Run a Python Script in ABAQUS:

### Within ABAQUS/CAE:

- Open ABAQUS/CAE.
  - Navigate to **File > Run Script.**
  - Browse to your .py file, select it, and click OK to run the script.
- Use the ABAQUS execution command if running the script outside CAE:
```console
 abaqus job=job_name user=subroutine_name script=script_name.py
```

