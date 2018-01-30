
#include "itkSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageMaskSpatialObject.h"

#include "itkMeanSquaresImageToImageMetricv4.h"

#include "itkLBFGSBOptimizerv4.h"
#include "itkImageRegistrationMethodv4.h"

#include "itkVersorRigid3DTransform.h"

#include "itkCommand.h"

#include "itkIdentityTransform.h"

int main( int argc, char *argv[] )
{

  const    unsigned int    ImageDimension = 3;
  typedef  double           PixelType;

  typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
  typedef itk::Image< PixelType, ImageDimension >  MovingImageType;
  typedef itk::Image< unsigned char, ImageDimension > MaskImageType;
  typedef itk::ImageMaskSpatialObject< ImageDimension >  FixedImageMaskType;
  typedef itk::IdentityTransform< double, ImageDimension > IdentityTransformType;

  typedef itk::MeanSquaresImageToImageMetricv4<
                                          FixedImageType,
                                          MovingImageType >    MetricType;


  typedef itk::LBFGSBOptimizerv4       OptimizerType;

  typedef itk::ImageRegistrationMethodv4<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  typedef itk::VersorRigid3DTransform< double > TransformType;

  TransformType::Pointer transform = TransformType::New();

  std::string outputFile(argv[1]);

  int numOfImages  = 0;
  std::sscanf(argv[2], "%d", &numOfImages);
  std::cout << "There are " << numOfImages << " images in the given series!" << std::endl;

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >   MaskImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  MaskImageReaderType::Pointer  maskImageReader  = MaskImageReaderType::New();

  FixedImageMaskType::Pointer spatialObjectMask = FixedImageMaskType::New();
  MetricType::Pointer         metric        = MetricType::New();


  maskImageReader->SetFileName( argv[3] );

  std::ofstream outfile_metric;
  outfile_metric.open(outputFile.c_str(), std::ofstream::out|std::ofstream::app);

  for (int imageIndex=1; imageIndex<numOfImages; imageIndex++) {

	  movingImageReader->SetFileName(  argv[4+imageIndex] );
	  fixedImageReader->SetFileName( argv[4] );

      OptimizerType::Pointer      optimizer     = OptimizerType::New();
      RegistrationType::Pointer   registration  = RegistrationType::New();

      registration->SetMetric(        metric        );
      registration->SetOptimizer(     optimizer     );
      registration->SetInitialTransform(transform);

	  FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
	  MovingImageType::ConstPointer movingImage = movingImageReader->GetOutput();
	  MaskImageType::ConstPointer maskImage = maskImageReader->GetOutput();

	  fixedImageReader->Update();
	  movingImageReader->Update();
	  maskImageReader->Update();

	  metric->SetMovingImage(movingImage);
	  metric->SetFixedImage(fixedImage);

	  std::cout << argv[4] << std::endl;
	  std::cout << argv[4+imageIndex] << std::endl;

	  spatialObjectMask->SetImage( maskImage );
	  //metric->SetMovingImageMask(spatialObjectMask);
	  //metric->SetFixedImageMask(spatialObjectMask);

	  optimizer->SetMaximumNumberOfCorrections(5);
	  optimizer->SetMaximumNumberOfFunctionEvaluations(5);

	  std::cout << metric->GetCurrentValue() << std::endl;
	  std::cout << metric->GetNumberOfValidPoints() << std::endl;


	  outfile_metric << metric->GetCurrentValue() << std::endl;


  }

  outfile_metric.close();
}
