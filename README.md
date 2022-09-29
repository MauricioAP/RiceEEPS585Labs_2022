# RiceEEPS585Labs
## Labs for the Computational Geophysical section of EEPS585
### Mauricio Araya, PhD.
### 09/06/2022

---

**Labs** setup given during class, and posted in canvas. For questions about software tools, please read file Tools.txt

---
**Lab 1:** Acquisition (09/16/2022)

1.1 create 2 different onshore-like acquistion plans (using acquisition.py)

1.2 run the illumination tool (illumination2D.py) with the SEG model for the 2 acquisions of 1.1

1.3 visualize the corresponding illumination maps (ill_view.py)

1.4 *Bonus* repeat experiment till the ideal acquisition to illuminate base of salt is achieve

The report needs to be delivered by **10/16/2022**, it should contain description of the results for item (1-3), *bonus* is optional and if achieved correctly it will be added to the final exam score.

---
**Lab 2:** Seismic Modeling and Imaging (09/30/2022)
Modeling, run the wave modeling application (modeling2D.py) on two different models (simple layered and SEG/EAGE).

2.1 For layered model (3 layers, 300x400) with a shot located in the middle of the surface,
    compare analytical and experimental travel time to selected reflectors.
    
2.2 For complex synth model (SEG/EAGE salt model) with shot from any location on the surface,
    take screen shots every n different times and make sense of it (explanation in terms of wavefront coherence, amplitude, location).

Imaging (RTM), run the migration tool (migration2D.py) on the models provided.

2.3 For layer model (3 layers, 100x200) check location of main reflectors.

2.4 For complex synth model (SEG/EAGE salt model) explain the image in terms of main reflectors and shadows (lack of reflectors) 
 in comparison to the velocity model.
 
2.5 Try different combinations of numbers of shots, iterations and peak frequency.

Bonus: if you create/use a different velocity models and repeat and report the experiments, this will grant 20% of the exam.

The report needs to be delivered by **11/04/2022**, it should contain description of the results for items (1-5), *bonus* is optional and if achieved correctly it will be added to the final exam score.
