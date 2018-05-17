This is a software built from ITK for temporal registration as described in http://people.csail.mit.edu/ruizhi/tmpReg. Temporal rigid registration is released. 

Usage:
${binPath}TemporalRegistration_Rigid outputPath numOfImages fixedImageMask weightOfTemporalRotationSmoothness weightOfTemporalTranslationSmoothness fixedImage movingImage1 [...]

Example:
${binPath}TemporalRegistration_Rigid ../registration_results/ 5 ../imageDir/1_mask.nii.gz 1 0 ../imageDir/1.nii.gz ../imageDir/2.nii.gz ../imageDir/3.nii.gz ../imageDir/4.nii.gz ../imageDir/5.nii.gz

Please contact ruizhi [at] mit.edu for any question.

-----
Links
-----

* Homepage: https://itk.org
* Download: https://itk.org/ITK/resources/software.html
* Mailing List: https://itk.org/ITK/help/mailing.html
* Book: https://itk.org/ITK/help/book.html
* Help: https://itk.org/ITK/help/help.html
* Examples: https://itk.org/ITKExamples/
* Bugtracker: https://issues.itk.org/
* Submit a patch: https://itk.org/Wiki/ITK/Git/Develop

-----
About ITK
-----

ITK is an open-source, cross-platform C++ toolkit for segmentation and
registration. Segmentation is the process of identifying and classifying
data found in a digitally sampled representation. Typically the sampled
representation is an image acquired from such medical instrumentation as
CT or MRI scanners. Registration is the task of aligning or developing
correspondences between data. For example, in the medical environment, a
CT scan may be aligned with a MRI scan in order to combine the information
contained in both.

The toolkit may be built from source using CMake (http://cmake.org).

