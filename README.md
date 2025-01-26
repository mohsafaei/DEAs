## In silico actuation performance investigation of dielectric elastomers with TPMS geometries [(Link)](https://doi.org/10.1016/j.euromechsol.2024.105540)


Part of the code for creating TPMS structures was developed by:<br/>
- **Fayyaz Nosouhi**: [email](dehnavifn@gmail.com), **Saeed Khaleghi**: [email](saeedkhaleghi123@gmail.com) <br/>

The remaining parts of the code, including the definition of DEA composites and their implementation in ABAQUS, were developed by:
- **Mohammad Ali Safaei** [**Email**](mohammadsf1998@gmail.com), [**GoogleScholar**](https://scholar.google.com/citations?user=jD_-4JcAAAAJ&hl=fa) , [**Linkedin**](https://www.linkedin.com/in/mohsafaei) <br />  

This Python script requires the use of a UEL subroutine, which implements the constitutive equations for a DEA element. In this regard, the UEL subroutine developed by Ehsan Hajiesmaili was utilized.    
> The UEL file and is available in the supplementary material of the following research article: [Link](https://pubs.aip.org/aip/jap/article/129/15/151102/1025587/Dielectric-elastomer-actuators)  
This project was done in **University of Tehran, December 2024**



## Steps to Run a Python Script in ABAQUS:

### Within ABAQUS/CAE:

- Open ABAQUS/CAE.
  - Navigate to **File > Run Script.**
  - Browse to your .py file, select it, and click OK to run the script.
- Use the ABAQUS execution command if running the script outside CAE:
```python
 abaqus job=job_name user=subroutine_name script=script_name.py
```

