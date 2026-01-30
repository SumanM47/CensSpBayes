Installing CensSpBayes
---------------------------------------------------------------------

# Step 1: Devtools

You would require the R package <tt>devtools</tt> to proceed with installation. If you do not have <tt>devtools</tt>, please install <tt>devtools</tt> first as below. If you have <tt>devtools</tt> installed in your machine, proceed to the next step.

<li><tt>install.packages("devtools")</tt></li>

# Step 2: INLA

This package requires the latest version of <tt>INLA</tt> package as a dependency. You can check out detailed installation instructions at <https://www.r-inla.org/download-install>. A brief description is given below. If you have the latest version of <tt>INLA</tt> on your machine, proceed to the next step.

If you are on Windows or Linux:
Use the following code to install <tt>INLA</tt> on your machine.

<li><tt>install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)</tt></li>


If you are on MacOS:
Please check the above website for the latest guidelines. Currently the following is required to install <tt>INLA</tt> on MacOS machines using R from [Homebrew](https://brew.sh/).
<li><tt>opt <- options()</tt></li>
<li><tt>options(pkgType="both")</tt></li>
<li><tt>install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)</tt></li>
<li><tt>options(opt)</tt></li>


# Step 3: Installing CensSpBayes 

With <tt>devtools</tt> and <tt>INLA</tt> on your machine, you can load it and install the <tt>CensSpBayes</tt> in the following way:

<li><tt>library(devtools)</tt></li>
<li><tt>devtools::install_github("SumanM47/CensSpBayes")</tt></li>

# Step 4: Using CensSpBayes

Please see the vignette at <https://sumanm47.github.io/Vignettes/CensSpBayes/>.
