# Genome-simulator
This is a web application using 7 models of GTR(general time reversible model) family to simulate genome. The models are GTR, K80, HKY, TN93, JC69, F81 and F84.

## Description
In EasyGTR, the web application will ask users to first choose the type of model that used to simulate Genome. The unique parameters of each model will show up automatically after selection. The application provides default values of parameters for each model. Besides these, users also need to set the insertion rate, deletion rate. length of root sequence and the number of generated random trees. EasyGTR also provides default values for all parameters.  
If not choosing by hand, users can also upload the model type and their parameters by uploading a text file. The format of the file is fixed. Users could find more details of the format in the example_input.txt file provided on github website.  

## Usage
The website of the web application is [easyGTR](https://xinning.shinyapps.io/easyGTR/)    
The original page of easyGTR should be like this:  
![image](https://github.com/luanxinning/Genome-simulator/assets/90717695/ef908f4f-e029-4ad7-a23a-ffb4eecfe1a0)  


  
If you want to upload the parameters from a file, choose "upload" in model type firstly. Click "browse" to upload the input file (check the format!) Then click "run".  
If you want to choose parameters by hand, firstly choose a model and change its parameters. Click "run" to get results.   
![image](https://github.com/luanxinning/Genome-simulator/assets/90717695/1fd3d8f8-5daa-4938-8482-e2bd7ed6069d)

  
 



## Results
The results parts are seperated into 6 tabs(Deletion, Insertion, Parameters, Site Rates Plot, Subtition and Simulation). Deletion tab includes a plot for deleltion size of distribution, and a detailed summary for deletion. Inerstion tab is similar with Deletion tab. In parameter tab ,this is a rate of matrix plot and a plot for equilibrium distribution. It also has a summary for rate parameters. Site Rates Plot tab has a plot for total site rate according to their positions; A phylogenetic tree is included in Subtition tab; In the last tab Simulation, there is a large plot will colored sites for each sequence. Meanwhile, the txt result of aligned sequences could be downloaded from Simulation tab.  
![image](https://github.com/luanxinning/Genome-simulator/assets/90717695/19966e79-ec87-4096-9da5-9bff93ee4361)

