# ElemeNT
ElemeNT is an interactive tool (implemented in Perl) for rapid and convenient detection of core promoter elements within a given sequence.
It can be executed via the following web interface: [ElemeNT 2023](https://www.juven-gershonlab.org/resources/element-v2023/) or executed via command line. 
The content of the folders located under ElemeNT2023 include:  
**ElemeNT/ElemeNT2023/installPackage** - zip of the installation package of ElemeNT2023 command line vesion.  
**ElemeNT/ElemeNT2023/src/ElemeNT** - the sources of ElemeNT2023.  
**ElemeNT/ElemeNT2023/src/Manuscript** - the source code used for the analysis included in the manuscript   
Below are the installation instruction for ElemeNT command line package: 
# Installation instructions
Install ElemeNT by extracting the ElemeNT_V2023.zip file into your chosen destination folder.  
## Content of the zip package: 
**Software executables:**  
ElemeNT_V2023.exe - for Windows OS (built under Windows 64bit).  
ElemeNT_V2023_binary - for Unix OS (built under Red Hat Enterprise Linux Server release 7.9 (Maipo))  
To allow execution, the file should be granted with execution permission.   
**The Elements folder**:  which contains the elements’ position weight matrices (PWM) files. The nucleotide order in the PWM files (top to bottom) is A C G T.  
Note – if needed, each PWM file’s contents can be modified to search for position weight matrices other than the defaults.  
**Config.txt file**- contains the parameter default values for the ElemeNT run. 
       The file name should not be changed. The file content is case sensitive.
Make sure you leave one blank space after each “:” character in the config.txt file. 
See “Configuration Settings” section for detailed explanations.  
**Sample_input.txt** - is an example of input file from which the sequence/s will be read   
**OutputElemeNT_V2023.txt** - is an example output of ElemeNT run using the sequences included within sample_input.txt as an input.   
**ElemeNT_v2023_UserGuide.pdf** - includes further details on the product and its parameter options.


