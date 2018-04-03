// Software Guide : BeginLatex
//
//
// Software Guide : EndLatex

#include "itkImageRegistrationMethodv4.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "itkBSplineTransform.h"
#include "itkVersorRigid3DTransform.h"

#include "itkLBFGSBOptimizerv4.h"

#include "itkSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkBSplineTransformInitializer.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkWarpImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "itkANTSNeighborhoodCorrelationImageToImageTemporalMetricv4.h"

#include "itkCommand.h"

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{

public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );

protected:
  RegistrationInterfaceCommand() {};

  // For convenience, we declare types useful for converting pointers
  // in the \code{Execute()} method.
  //
public:
  typedef   TRegistration      RegistrationType;
  typedef   RegistrationType * RegistrationPointer;
  typedef itk::LBFGSBOptimizerv4       OptimizerType;
  typedef   OptimizerType * OptimizerPointer;

  // Two arguments are passed to the \code{Execute()} method: the first
  // is the pointer to the object which invoked the event and the
  // second is the event that was invoked.
  //
  void Execute( itk::Object * object,
                const itk::EventObject & event) ITK_OVERRIDE
    {
    // First we verify that the event invoked is of the right type,
    // \code{itk::MultiResolutionIterationEvent()}.
    // If not, we return without any further action.
    //
    if( !(itk::MultiResolutionIterationEvent().CheckEvent( &event ) ) )
      {
      return;
      }
    // We then convert the input object pointer to a RegistrationPointer.
    // Note that no error checking is done here to verify the
    // \code{dynamic\_cast} was successful since we know the actual object
    // is a registration method. Then we ask for the optimizer object
    // from the registration method.
    //
    RegistrationPointer registration = static_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer =  static_cast< OptimizerPointer >( registration->GetModifiableOptimizer() );

    unsigned int currentLevel = registration->GetCurrentLevel();
    typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactors = registration->GetShrinkFactorsPerDimension( currentLevel );
    typename RegistrationType::SmoothingSigmasArrayType smoothingSigmas = registration->GetSmoothingSigmasPerLevel();

    std::cout << "-------------------------------------" << std::endl;
    std::cout << " Current level = " << currentLevel << std::endl;
    std::cout << "    shrink factor = " << shrinkFactors << std::endl;
    std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    try {
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << std::endl;
    }
    catch(itk::ExceptionObject)
    {
    }

    if ( registration->GetCurrentLevel() == 0 )
      {
      }
    else
      {
    	optimizer->SetInitialPosition(optimizer->GetCurrentPosition());
      }
    }

  // Another version of the \code{Execute()} method accepting a \code{const}
  // input object is also required since this method is defined as pure virtual
  // in the base class.  This version simply returns without taking any action.
  //
  void Execute(const itk::Object * , const itk::EventObject & ) ITK_OVERRIDE
    {
    return;
    }
};
//  The following section of code implements an observer
//  that will monitor the evolution of the registration process.
//
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate(): m_CumulativeIterationIndex(0) {};

public:
  typedef itk::LBFGSBOptimizerv4       OptimizerType;
  typedef   const OptimizerType *                               OptimizerPointer;
  std::string outputOptimizationLog;

  void SetOutputOptimizationLog (std::string outputFile)
  {
	  this->outputOptimizationLog = outputFile;
  }

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
  {
	  Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
  {
	  OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
	  if( !(itk::IterationEvent().CheckEvent( &event )) )
	  {
		  return;
	  }
	  std::cout << optimizer->GetCurrentIteration() << "   ";
	  std::cout << optimizer->GetValue() << "   ";
	  std::cout << m_CumulativeIterationIndex++ << std::endl;
  }
private:
  unsigned int m_CumulativeIterationIndex;
};

int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
	  std::cerr << "Missing Parameters " << std::endl;
	  std::cerr << "Usage: " << argv[0];
	  std::cerr << " fixedImageFile  movingImageFile outputImagefile  ";
	  std::cerr << " [differenceOutputfile] [differenceBeforeRegistration] ";
	  std::cerr << " [deformationField] ";
	  return EXIT_FAILURE;
  }

  const    unsigned int    ImageDimension = 3;
  typedef  double           PixelType;

  typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
  typedef itk::Image< PixelType, ImageDimension >  MovingImageType;
  typedef itk::Image< unsigned char, ImageDimension > MaskImageType;
  typedef itk::ImageMaskSpatialObject< ImageDimension >  FixedImageMaskType;

  typedef double CoordinateRepType;

  const unsigned int SpaceDimension = ImageDimension;
  const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;

  typedef itk::BSplineTransform<
                            CoordinateRepType,
                            SpaceDimension,
                            SplineOrder >     TransformType;

  typedef itk::LBFGSBOptimizerv4       OptimizerType;

  typedef itk::ANTSNeighborhoodCorrelationImageToImageTemporalMetricv4<
                                          FixedImageType,
                                          MovingImageType >    MetricType;

  typedef itk::ImageRegistrationMethodv4<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  std::string outputFolder(argv[1]);

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

  OptimizerType::ParametersType transformParameters;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );

  maskImageReader->SetFileName( argv[3] );

  std::ofstream outfile_metric;
  std::string outputMetricValues = outputFolder + "metric_values.txt";
  outfile_metric.open(outputMetricValues.c_str(), std::ofstream::out|std::ofstream::app);

  std::ofstream outfile_parameters;
  std::string outputTransformParameters = outputFolder + "transform_parameters.txt";
  outfile_parameters.open(outputTransformParameters.c_str(), std::ofstream::out|std::ofstream::app);

  double w1;
  std::sscanf(argv[4], "%lf", &w1);
  double w2;
  std::sscanf(argv[5], "%lf", &w2);
  double zero=0;

  std::string transformModel = "bspline";
  metric->SetTransformModel(transformModel);

  TransformType::Pointer  outputTransform = TransformType::New();
  TransformType::Pointer  initialTransform = TransformType::New();

  for (int imageIndex=1; imageIndex<numOfImages; imageIndex++) {

	  movingImageReader->SetFileName(  argv[6] );
	  fixedImageReader->SetFileName( argv[6+imageIndex] );

	  std::string movingImageName(argv[6]);
	  std::string slash = "/";
	  std::size_t movingSlashIndex = movingImageName.find_last_of(slash);
	  std::string fixedImageName(argv[6+imageIndex]);
	  std::size_t fixedSlashIndex = fixedImageName.find_last_of(slash);

	  std::string movedImageName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + ".nii.gz";
	  std::cout << "Moved Image Name: " << movedImageName << std::endl;
	  std::string maskedMovedImageName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_masked.nii.gz";
	  std::cout << "Moved Image Name: " << maskedMovedImageName << std::endl;
	  std::string warpFieldName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_warp.nii.gz";
	  std::cout << "Warp Image Name: " << warpFieldName << std::endl;

	  std::string inverseMovedImageName = outputFolder + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_to_" + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + ".nii.gz";
	  std::cout << "Inverse Moved Image Name: " << inverseMovedImageName << std::endl;
	  std::string inverseWarpFieldName = outputFolder + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_to_" + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_inverseWarp.nii.gz";
	  std::cout << "Inverse Warp Image Name: " << inverseWarpFieldName << std::endl;


	  FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
	  MovingImageType::ConstPointer movingImage = movingImageReader->GetOutput();
	  MaskImageType::ConstPointer maskImage = maskImageReader->GetOutput();

	  fixedImageReader->Update();
	  movingImageReader->Update();
	  maskImageReader->Update();

	  registration->SetFixedImage(  fixedImage   );
	  registration->SetMovingImage(   movingImage   );
	  spatialObjectMask->SetImage( maskImage );
	  metric->SetMovingImageMask(spatialObjectMask);

	  outfile_parameters << outputTransform->GetParameters() << std::endl;

	  // Initialize the transform
	  typedef itk::BSplineTransformInitializer< TransformType,
												FixedImageType>      InitializerType;
	  InitializerType::Pointer transformInitializer = InitializerType::New();

	  TransformType::MeshSizeType             meshSize;
	  meshSize.SetElement(0, 12);
	  meshSize.SetElement(1, 12);
	  meshSize.SetElement(2, 9);

	  if (imageIndex==1) {
		  metric->SetTemporalSmoothness1(zero);
		  metric->SetTemporalSmoothness2(zero);

		  transformInitializer->SetTransform(   initialTransform );
		  transformInitializer->SetImage(  fixedImageReader->GetOutput() );
		  transformInitializer->SetTransformDomainMeshSize( meshSize );
		  transformInitializer->InitializeTransform();

		  registration->SetInitialTransform( initialTransform );
		  registration->InPlaceOn();
	  } else {
		  metric->SetTemporalSmoothness1(w1);
		  metric->SetTemporalSmoothness2(zero);

		  outfile_parameters << outputTransform->GetParameters() << std::endl;

		  transformInitializer->SetTransform(   outputTransform );
		  transformInitializer->SetImage(  fixedImageReader->GetOutput() );
		  transformInitializer->SetTransformDomainMeshSize( meshSize );
		  transformInitializer->InitializeTransform();

		  registration->SetInitialTransform( outputTransform );
		  registration->InPlaceOn();
	  }

	  const unsigned int numParameters = initialTransform->GetNumberOfParameters();
	  double* previousTransformParameters = new double[numParameters];
	  if (imageIndex==1) {
		  for (int i=0; i<int(numParameters); i++)
		  {
			  previousTransformParameters[i] = 0;
		  }
		  metric->SetPreviousTransformParameters(previousTransformParameters, numParameters);
	  }

	  const unsigned int numberOfLevels3 = 3;

	  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	  shrinkFactorsPerLevel.SetSize( numberOfLevels3 );
	  shrinkFactorsPerLevel[0] = 4;
	  shrinkFactorsPerLevel[1] = 2;
	  shrinkFactorsPerLevel[2] = 1;
	  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	  smoothingSigmasPerLevel.SetSize( numberOfLevels3 );
	  smoothingSigmasPerLevel[0] = 4;
	  smoothingSigmasPerLevel[1] = 2;
	  smoothingSigmasPerLevel[2] = 0;

	  registration->SetNumberOfLevels( numberOfLevels3 );
	  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
	  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

	  //  Next we set the parameters of the LBFGSB Optimizer. Note that
	  //  this optimizer does not support scales estimator and sets all
	  //  the parameters scales to one.
	  //  Also, we should set the boundary condition for each variable, where
	  //  \code{boundSelect[i]} can be set as: \code{UNBOUNDED},
	  //  \code{LOWERBOUNDED}, \code{BOTHBOUNDED}, \code{UPPERBOUNDED}.
	  //

	  std::cout << "Num of transform parameters:" << numParameters << std::endl;

	  OptimizerType::BoundSelectionType boundSelect( numParameters );
	  OptimizerType::BoundValueType upperBound( numParameters );
	  OptimizerType::BoundValueType lowerBound( numParameters );

	  boundSelect.Fill( OptimizerType::UNBOUNDED );
	  upperBound.Fill( 0.0 );
	  lowerBound.Fill( 0.0 );

	  optimizer->SetBoundSelection( boundSelect );
	  optimizer->SetUpperBound( upperBound );
	  optimizer->SetLowerBound( lowerBound );
	  optimizer->SetCostFunctionConvergenceFactor( 1e+12 );
	  optimizer->SetGradientConvergenceTolerance( 1.0e-35 );
	  optimizer->SetNumberOfIterations( 30 );
	  optimizer->SetMaximumNumberOfFunctionEvaluations( 100 );
	  optimizer->SetMaximumNumberOfCorrections( 5 );
	  if (imageIndex==1) {
	  	  optimizer->SetInitialPosition(initialTransform->GetParameters());
	  } else {
	   	  optimizer->SetInitialPosition(outputTransform->GetParameters());
	  }

	  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	  optimizer->AddObserver( itk::IterationEvent(), observer );

	  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
	  CommandType::Pointer command = CommandType::New();
	  registration->AddObserver( itk::MultiResolutionIterationEvent(), command );

	  // Add time and memory probes
	  itk::TimeProbesCollectorBase chronometer;
	  itk::MemoryProbesCollectorBase memorymeter;

	  std::cout << std::endl << "Starting Registration" << std::endl;

	  try
		{
		memorymeter.Start( "Registration" );
		chronometer.Start( "Registration" );

		registration->Update();

		chronometer.Stop( "Registration" );
		memorymeter.Stop( "Registration" );

		optimizer->RemoveAllObservers();
		registration->RemoveAllObservers();

		std::cout << "Optimizer stop condition = "
				  << registration->GetOptimizer()->GetStopConditionDescription()
				  << std::endl;
		}
	  catch( itk::ExceptionObject & err )
		{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
		}

	  double metricValue = optimizer->GetCurrentMetricValue();
	  outfile_metric << metricValue <<std::endl;

	  if (imageIndex==1) {
		  outputTransform = initialTransform->Clone();
	  }
	  transformParameters = outputTransform->GetParameters();

	  for (int i=0; i<int(numParameters); i++)
	  {
		  previousTransformParameters[i] = transformParameters[i];
	  }
	  metric->SetPreviousTransformParameters(previousTransformParameters, numParameters);

	  outfile_parameters << movedImageName << std::endl;
	  outfile_parameters << transformParameters << std::endl;

	  // Finally we use the last transform in order to resample the image.
	  //
	  typedef itk::ResampleImageFilter<
								MovingImageType,
								FixedImageType >    ResampleFilterType;

	  ResampleFilterType::Pointer resample = ResampleFilterType::New();

	  resample->SetTransform( outputTransform );
	  resample->SetInput( movingImageReader->GetOutput() );

	  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
	  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
	  resample->SetOutputSpacing( fixedImage->GetSpacing() );
	  resample->SetOutputDirection( fixedImage->GetDirection() );
	  resample->SetDefaultPixelValue( 0 );

	  //typedef  unsigned char  OutputPixelType;
	  typedef PixelType OutputPixelType;
	  typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;

	  typedef itk::CastImageFilter<
							FixedImageType,
							OutputImageType > CastFilterType;

	  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	  WriterType::Pointer      writer =  WriterType::New();
	  CastFilterType::Pointer  caster =  CastFilterType::New();

	  writer->SetFileName( movedImageName );

	  caster->SetInput( resample->GetOutput() );
	  writer->SetInput( caster->GetOutput()   );

	  try
		{
		writer->Update();
		}
	  catch( itk::ExceptionObject & err )
		{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
		}

	  // Generate the explicit deformation field resulting from
	  // the registration.
	  typedef itk::Vector< float, ImageDimension >          VectorPixelType;
	  typedef itk::Image< VectorPixelType, ImageDimension > DisplacementFieldImageType;

	  typedef itk::TransformToDisplacementFieldFilter<
							DisplacementFieldImageType,
							CoordinateRepType >             DisplacementFieldGeneratorType;

	  /** Create an setup displacement field generator. */
	  DisplacementFieldGeneratorType::Pointer dispfieldGenerator =
													 DisplacementFieldGeneratorType::New();

	  DisplacementFieldImageType::Pointer resampledInverseDispField = DisplacementFieldImageType::New();

	  dispfieldGenerator->UseReferenceImageOn();
	  dispfieldGenerator->SetReferenceImage( fixedImage );
	  dispfieldGenerator->SetTransform( outputTransform );
	  try
		{
		dispfieldGenerator->Update();
		}
	  catch ( itk::ExceptionObject & err )
		{
		std::cerr << "Exception detected while generating deformation field";
		std::cerr << " : "  << err << std::endl;
		return EXIT_FAILURE;
		}

	  typedef itk::ImageFileWriter< DisplacementFieldImageType >  FieldWriterType;
	  FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

	  DisplacementFieldImageType::Pointer dispField = dispfieldGenerator->GetOutput();

	  DisplacementFieldImageType::RegionType dispFieldRegion = dispField->GetLargestPossibleRegion();
	  itk::ImageRegionIteratorWithIndex< DisplacementFieldImageType > it( dispField, dispFieldRegion );

	  VectorPixelType zeroVector;
	  for (int i=0; i<(int)ImageDimension; i++) {
		  zeroVector[i] = 0.0;
	  }
	  it.GoToBegin();
	  while( !it.IsAtEnd() ) {
		if(maskImage->GetPixel(it.GetIndex()) == 0) {
			it.Set(zeroVector);
		}
	    ++it;
	  }

	  fieldWriter->SetInput( dispField );
	  fieldWriter->SetFileName( warpFieldName );
	  try
	    {
		fieldWriter->Update();
	    }
	  catch( itk::ExceptionObject & excp )
	    {
		std::cerr << "Exception thrown " << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	    }

	  typedef itk::WarpImageFilter<
			  	  	  FixedImageType,
					  MovingImageType,
					  DisplacementFieldImageType >             WarpImageFilterType;

	  WarpImageFilterType::Pointer warpMovingImageFilter = WarpImageFilterType::New();

	  warpMovingImageFilter->SetInput(movingImage);
	  warpMovingImageFilter->SetOutputParametersFromImage(fixedImage);
	  warpMovingImageFilter->SetDisplacementField(dispField);
	  warpMovingImageFilter->Update();

	  WriterType::Pointer      writer2 =  WriterType::New();

	  writer2->SetFileName( maskedMovedImageName );
	  writer2->SetInput( warpMovingImageFilter->GetOutput() );

	  try
		{
		writer2->Update();
		}
	  catch( itk::ExceptionObject & err )
		{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
		}

	  // Finally we use the last inverse transform in order to resample the image.
	  //
	  DisplacementFieldImageType::Pointer inverseDispField = dispField;

	  DisplacementFieldImageType::RegionType inverseDispFieldRegion = inverseDispField->GetLargestPossibleRegion();
	  itk::ImageRegionIteratorWithIndex< DisplacementFieldImageType > itInverseDisp( inverseDispField, inverseDispFieldRegion );

	  VectorPixelType vector;
	  itInverseDisp.GoToBegin();
	  while( !itInverseDisp.IsAtEnd() )
	    {
		vector = itInverseDisp.Get();
		for (int i=0; i<(int)ImageDimension; i++) {
			vector[i] = vector[i] * -1.0;
		}
		itInverseDisp.Set( vector );
	    ++itInverseDisp;
	    }

	  typedef itk::ResampleImageFilter<
			  	  	  	  DisplacementFieldImageType,
						  DisplacementFieldImageType >    ResampleDispFieldFilterType;

	  ResampleDispFieldFilterType::Pointer dispFieldResample = ResampleDispFieldFilterType::New();
	  dispFieldResample->SetTransform( outputTransform );
	  dispFieldResample->SetInput( inverseDispField );
	  dispFieldResample->SetSize(    dispField->GetLargestPossibleRegion().GetSize() );
	  dispFieldResample->SetOutputOrigin(  dispField->GetOrigin() );
	  dispFieldResample->SetOutputSpacing( dispField->GetSpacing() );
	  dispFieldResample->SetOutputDirection( dispField->GetDirection() );
	  dispFieldResample->SetDefaultPixelValue( zeroVector );
	  typedef itk::CastImageFilter<
			  	  	  	  DisplacementFieldImageType,
					  	  DisplacementFieldImageType > CastDispFieldFilterType;
	  CastDispFieldFilterType::Pointer dispFieldCaster = CastDispFieldFilterType::New();
	  dispFieldCaster->SetInput(dispFieldResample->GetOutput());
	  dispFieldCaster->Update();
	  resampledInverseDispField = dispFieldCaster->GetOutput();

	  WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();

	  warpImageFilter->SetInput(fixedImage);
	  warpImageFilter->SetOutputParametersFromImage(movingImage);
	  warpImageFilter->SetDisplacementField(resampledInverseDispField);
	  warpImageFilter->Update();

	  WriterType::Pointer      writer3 =  WriterType::New();

	  writer3->SetFileName( inverseMovedImageName );
	  writer3->SetInput( warpImageFilter->GetOutput() );

	  try
		{
		writer3->Update();
		}
	  catch( itk::ExceptionObject & err )
		{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
		}

	  // Generate the explicit inverse deformation field resulting from
	  // the registration.
	  typedef itk::ImageFileWriter< DisplacementFieldImageType >  FieldWriterType;
	  FieldWriterType::Pointer fieldWriter2 = FieldWriterType::New();

	  fieldWriter2->SetInput( resampledInverseDispField );
	  fieldWriter2->SetFileName( inverseWarpFieldName );
	  try
	    {
		  fieldWriter2->Update();
	    }
	  catch( itk::ExceptionObject & excp )
	    {
		std::cerr << "Exception thrown " << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	    }

  }

  outfile_metric.close();
  outfile_parameters.close();

  return EXIT_SUCCESS;
}
