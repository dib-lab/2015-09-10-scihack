## Scihack!

Code from the DIB lab's first attempt at a science hack session. For now it has a pipeline to retrieve the
sample data; it also has the code that generated the random samples, though not the raw data they were
sampled from. To get up and running, clone the repo:

    git clone git@github.com:dib-lab/2015-09-10-scihack.git

This pipeline depends on [pydoit](http://pydoit.org/) and [pandas](http://pandas.pydata.org/).
I recommend anaconda for managing your python environment, which will include pandas already.
To get pydoit:

    pip install doit

Then, to get the data:

	./pipeline run
