// Software Guide : BeginLatex
//
//
//
// Software Guide : EndLatex

#include "itkImageRegistrationMethodv4.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkLBFGSBOptimizerv4.h"

#include "itkSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMaskImageFilter.h"
#include "itkImageMomentsCalculator.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkTransformToDisplacementFieldFilter.h"

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

public:
  typedef   TRegistration      RegistrationType;
  typedef   RegistrationType * RegistrationPointer;
  typedef 	itk::LBFGSBOptimizerv4       OptimizerType;
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
    OptimizerPointer optimizer =  static_cast< OptimizerPointer >(registration->GetModifiableOptimizer() );

    unsigned int currentLevel = registration->GetCurrentLevel();
    typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactors =
      registration->GetShrinkFactorsPerDimension( currentLevel );
    typename RegistrationType::SmoothingSigmasArrayType smoothingSigmas =
      registration->GetSmoothingSigmasPerLevel();

    std::cout << "-------------------------------------" << std::endl;
    std::cout << " Current level = " << currentLevel << std::endl;
    std::cout << "    shrink factor = " << shrinkFactors << std::endl;
    std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    try {
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
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
  typedef  itk::LBFGSBOptimizerv4       OptimizerType;
  typedef  const OptimizerType * OptimizerPointer;
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
	  std::cout << optimizer->GetCurrentPosition() << "   ";
	  std::cout << optimizer->GetCachedDerivative() << "   ";
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
  typedef  double          PixelType;

  typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
  typedef itk::Image< PixelType, ImageDimension >  MovingImageType;
  typedef itk::Image< unsigned char, ImageDimension > MaskImageType;
  typedef itk::ImageMaskSpatialObject< ImageDimension >  ImageMaskType;

  typedef double CoordinateRepType;
  typedef itk::VersorRigid3DTransform< double > TransformType;

  typedef itk::LBFGSBOptimizerv4       OptimizerType;

  typedef itk::ANTSNeighborhoodCorrelationImageToImageTemporalMetricv4<
                                          FixedImageType,
                                          MovingImageType >    MetricType;

  typedef itk::ImageRegistrationMethodv4<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  typedef itk::MaskImageFilter< MovingImageType, MaskImageType > MaskFileterType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >   MaskImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  MaskImageReaderType::Pointer  fixedMaskImageReader  = MaskImageReaderType::New();

  ImageMaskType::Pointer spatialObjectFixedMask = ImageMaskType::New();

  OptimizerType::ParametersType transformParameters;

  MaskFileterType::Pointer maskFilter = MaskFileterType::New();

  std::string outputFolder(argv[1]);

  int numOfImages = 0;
  std::sscanf(argv[2], "%d", &numOfImages);
  std::cout << "There are " << numOfImages << " images in the given series!" << std::endl;

  fixedMaskImageReader->SetFileName( argv[3] );

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

  double* t = new double[6];
  for (int i=0; i<6; i++)
  {
	  t[i] = 0;
  }
  metric->SetPreviousTransformParameters(t, 6);
  std::string transformModel = "rigid";
  metric->SetTransformModel(transformModel);

  TransformType::Pointer  outputTransform1 = TransformType::New();
  TransformType::Pointer  outputTransform2 = TransformType::New();
  TransformType::Pointer  outputTransform3 = TransformType::New();
  TransformType::Pointer  outputTransform = TransformType::New();
  TransformType::Pointer  inverseOutputTransform = TransformType::New();

  TransformType::Pointer  initialTransform1 = TransformType::New();
  TransformType::Pointer  initialTransform2 = TransformType::New();
  TransformType::Pointer  initialTransform3 = TransformType::New();

  for (int imageIndex=1; imageIndex<numOfImages; imageIndex++) {


	  movingImageReader->SetFileName(  argv[6+imageIndex] );
	  fixedImageReader->SetFileName( argv[6] );

	  std::string movingImageName(argv[6+imageIndex]);
	  std::string slash = "/";
	  std::size_t movingSlashIndex = movingImageName.find_last_of(slash);
	  std::string fixedImageName(argv[6]);
	  std::size_t fixedSlashIndex = fixedImageName.find_last_of(slash);

	  std::string movedImageName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + ".nii.gz";
	  std::cout << "Moved Image Name: " << movedImageName << std::endl;
	  std::string warpFieldName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_warp.nii.gz";
	  std::cout << "Warp Image Name: " << warpFieldName << std::endl;

	  std::string inverseMovedImageName = outputFolder + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_to_" + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + ".nii.gz";
	  std::cout << "Inverse Moved Image Name: " << inverseMovedImageName << std::endl;
	  std::string inverseWarpFieldName = outputFolder + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_to_" + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_inverseWarp.nii.gz";
	  std::cout << "Inverse Warp Image Name: " << inverseWarpFieldName << std::endl;

	  FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
	  MovingImageType::ConstPointer movingImage = movingImageReader->GetOutput();
	  MaskImageType::ConstPointer fixedMaskImage = fixedMaskImageReader->GetOutput();

	  fixedMaskImageReader->Update();
	  fixedImageReader->Update();
	  movingImageReader->Update();

	  registration->SetFixedImage(  fixedImage   );
	  registration->SetMovingImage(   movingImage   );
	  spatialObjectFixedMask->SetImage( fixedMaskImage );
	  metric->SetFixedImageMask(spatialObjectFixedMask);

	  maskFilter->SetInput( fixedImage );
	  maskFilter->SetMaskImage( fixedMaskImage );
	  MovingImageType::ConstPointer maskedFixedImage = maskFilter->GetOutput();
	  maskFilter->Update();

	  typedef itk::ImageMomentsCalculator< FixedImageType > FixedImageCalculatorType;
	  FixedImageCalculatorType::Pointer FixedImageCalculator = FixedImageCalculatorType::New();
	  FixedImageCalculator->SetImage( maskedFixedImage );
	  FixedImageCalculator->Compute();

	  FixedImageCalculatorType::VectorType movingCenter = FixedImageCalculator->GetCenterOfGravity();
	  std::cout << "Moving Center: " << movingCenter << std::endl;

	  outfile_parameters << outputTransform->GetParameters() << std::endl;

	  // Initialize the transform
	  typedef itk::CenteredTransformInitializer<
			  TransformType,
			  FixedImageType,
		      MovingImageType >  TransformInitializerType;
	  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
	  TransformInitializerType::Pointer initializer2 = TransformInitializerType::New();


	  // Three-level registrtaion
	  //
	  if (imageIndex==1) {
		  metric->SetTemporalSmoothness1(zero);
		  metric->SetTemporalSmoothness2(zero);

		  initializer->SetTransform(   initialTransform3 );
	      initializer->SetFixedImage(  fixedImage );
	      initializer->SetMovingImage( movingImage );
	      initializer->MomentsOn();
		  initializer->InitializeTransform();

		  typedef TransformType::VersorType  VersorType;
		  typedef VersorType::VectorType     VectorType;
		  VersorType     rotation;
		  VectorType     axis;
		  axis[0] = 0;
		  axis[1] = 0;
		  axis[2] = 1;
		  const double angle = 0;
		  rotation.Set(  axis, angle  );
		  initialTransform3->SetRotation( rotation );
		  initialTransform3->SetCenter(movingCenter);

		  registration->SetInitialTransform( initialTransform3 );
		  registration->InPlaceOn();
	  } else {
		  metric->SetTemporalSmoothness1(w1);
		  metric->SetTemporalSmoothness2(w2);

		  outputTransform3 = outputTransform->Clone();

		  initializer2->SetTransform(   outputTransform3 );
		  initializer2->SetFixedImage(  fixedImage );
		  initializer2->SetMovingImage( movingImage );
		  initializer2->InitializeTransform();

		  registration->SetInitialTransform( outputTransform3 );
		  registration->InPlaceOn();

		  std::cout << "outputTransform: " << outputTransform->GetParameters() << std::endl;
		  std::cout << "offset: " << outputTransform3->GetOffset() << std::endl;
		  std::cout << "center: " << outputTransform3->GetCenter() << std::endl;
		  std::cout << "outputTransform3: " << outputTransform3->GetParameters() << std::endl;
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
	  const unsigned int numParameters = outputTransform->GetNumberOfParameters();

	  OptimizerType::BoundSelectionType boundSelect( numParameters );
	  OptimizerType::BoundValueType upperBound( numParameters );
	  OptimizerType::BoundValueType lowerBound( numParameters );

	  boundSelect.Fill( OptimizerType::UNBOUNDED );
	  upperBound.Fill( 0.0 );
	  lowerBound.Fill( 0.0 );

	  optimizer->SetBoundSelection( boundSelect );
	  optimizer->SetUpperBound( upperBound );
	  optimizer->SetLowerBound( lowerBound );
	  optimizer->SetCostFunctionConvergenceFactor( 1e+7 );
	  optimizer->SetGradientConvergenceTolerance( 1.0e-35 );
	  optimizer->SetNumberOfIterations( 30 );
	  optimizer->SetMaximumNumberOfFunctionEvaluations( 200 );
	  optimizer->SetMaximumNumberOfCorrections( 6 );

	  if (imageIndex==1) {
	  	  optimizer->SetInitialPosition(initialTransform3->GetParameters());
	  } else {
	   	  optimizer->SetInitialPosition(transformParameters);
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

	  double metricValue3 = optimizer->GetCurrentMetricValue();
	  outfile_metric << metricValue3 << std::endl;

	  std::cout << initialTransform3->GetParameters() << std::endl;
	  std::cout << initialTransform3->GetMatrix() << std::endl;
	  std::cout << initialTransform3->GetTranslation() << std::endl;
	  std::cout << initialTransform3->GetRotationMatrix() << std::endl;

	  // Two-level registrtaion
	  //
	  if (imageIndex==1) {
		  metric->SetTemporalSmoothness1(zero);
		  metric->SetTemporalSmoothness2(zero);

		  initializer->SetTransform(   initialTransform2 );
	      initializer->SetFixedImage(  fixedImage );
	      initializer->SetMovingImage( movingImage );
	      initializer->MomentsOn();
		  initializer->InitializeTransform();

		  typedef TransformType::VersorType  VersorType;
		  typedef VersorType::VectorType     VectorType;
		  VersorType     rotation;
		  VectorType     axis;
		  axis[0] = 0;
		  axis[1] = 0;
		  axis[2] = 1;
		  const double angle = 0;
		  rotation.Set(  axis, angle  );
		  initialTransform2->SetRotation( rotation );
		  initialTransform2->SetCenter(movingCenter);

		  registration->SetInitialTransform( initialTransform2 );
		  registration->InPlaceOn();
	  } else {
		  metric->SetTemporalSmoothness1(w1);
		  metric->SetTemporalSmoothness2(w2);

		  outputTransform2 = outputTransform->Clone();

		  initializer->SetTransform(   outputTransform2 );
	      initializer->SetFixedImage(  fixedImage );
	      initializer->SetMovingImage( movingImage );
	      initializer->GeometryOn();
		  initializer->InitializeTransform();

		  registration->SetInitialTransform( outputTransform2 );
		  registration->InPlaceOn();

		  std::cout << "outputTransform2: " << outputTransform2->GetParameters() << std::endl;

	  }

	  const unsigned int numberOfLevels2 = 2;

	  shrinkFactorsPerLevel.SetSize( numberOfLevels2 );
	  shrinkFactorsPerLevel[0] = 2;
	  shrinkFactorsPerLevel[1] = 1;
	  smoothingSigmasPerLevel.SetSize( numberOfLevels2 );
	  smoothingSigmasPerLevel[0] = 2;
	  smoothingSigmasPerLevel[1] = 0;

	  registration->SetNumberOfLevels( numberOfLevels2 );
	  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
	  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

	  //  Next we set the parameters of the LBFGSB Optimizer. Note that
	  //  this optimizer does not support scales estimator and sets all
	  //  the parameters scales to one.
	  //  Also, we should set the boundary condition for each variable, where
	  //  \code{boundSelect[i]} can be set as: \code{UNBOUNDED},
	  //  \code{LOWERBOUNDED}, \code{BOTHBOUNDED}, \code{UPPERBOUNDED}.
	  //
	  boundSelect.Fill( OptimizerType::UNBOUNDED );
	  upperBound.Fill( 0.0 );
	  lowerBound.Fill( 0.0 );

	  optimizer->SetBoundSelection( boundSelect );
	  optimizer->SetUpperBound( upperBound );
	  optimizer->SetLowerBound( lowerBound );
	  optimizer->SetCostFunctionConvergenceFactor( 1e+7 );
	  optimizer->SetGradientConvergenceTolerance( 1.0e-35 );
	  optimizer->SetNumberOfIterations( 30 );
	  optimizer->SetMaximumNumberOfFunctionEvaluations( 200 );
	  optimizer->SetMaximumNumberOfCorrections( 6 );

	  if (imageIndex==1) {
	  	  optimizer->SetInitialPosition(initialTransform2->GetParameters());
	  } else {
	   	  optimizer->SetInitialPosition(transformParameters);
	  }

	  observer = CommandIterationUpdate::New();
	  optimizer->AddObserver( itk::IterationEvent(), observer );

	  command = CommandType::New();
	  registration->AddObserver( itk::MultiResolutionIterationEvent(), command );

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

	  double metricValue2 = optimizer->GetCurrentMetricValue();
	  outfile_metric << metricValue2 << std::endl;

	  // Single-level registration
	  //
	  if (imageIndex==1) {
		  metric->SetTemporalSmoothness1(zero);
		  metric->SetTemporalSmoothness2(zero);

		  initializer->SetTransform(   initialTransform1 );
	      initializer->SetFixedImage(  fixedImage );
	      initializer->SetMovingImage( movingImage );
	      initializer->MomentsOn();
		  initializer->InitializeTransform();

		  typedef TransformType::VersorType  VersorType;
		  typedef VersorType::VectorType     VectorType;
		  VersorType     rotation;
		  VectorType     axis;
		  axis[0] = 0;
		  axis[1] = 0;
		  axis[2] = 1;
		  const double angle = 0;
		  rotation.Set(  axis, angle  );
		  initialTransform1->SetRotation( rotation );
		  initialTransform1->SetCenter(movingCenter);

		  registration->SetInitialTransform( initialTransform1 );
		  registration->InPlaceOn();
	  } else {
		  metric->SetTemporalSmoothness1(w1);
		  metric->SetTemporalSmoothness2(w2);

		  outputTransform1 = outputTransform->Clone();

		  initializer->SetTransform(   outputTransform1 );
	      initializer->SetFixedImage(  fixedImage );
	      initializer->SetMovingImage( movingImage );
	      initializer->GeometryOn();
		  initializer->InitializeTransform();

		  registration->SetInitialTransform( outputTransform1 );
		  registration->InPlaceOn();

		  std::cout << "outputTransform1: " << outputTransform1->GetParameters() << std::endl;

	  }

	  const unsigned int numberOfLevels1 = 1;

	  shrinkFactorsPerLevel.SetSize( numberOfLevels1 );
	  shrinkFactorsPerLevel[0] = 1;
	  smoothingSigmasPerLevel.SetSize( numberOfLevels1 );
	  smoothingSigmasPerLevel[0] = 0;

	  registration->SetNumberOfLevels( numberOfLevels1 );
	  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
	  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

	  //  Next we set the parameters of the LBFGSB Optimizer. Note that
	  //  this optimizer does not support scales estimator and sets all
	  //  the parameters scales to one.
	  //  Also, we should set the boundary condition for each variable, where
	  //  \code{boundSelect[i]} can be set as: \code{UNBOUNDED},
	  //  \code{LOWERBOUNDED}, \code{BOTHBOUNDED}, \code{UPPERBOUNDED}.
	  //
	  boundSelect.Fill( OptimizerType::UNBOUNDED );
	  upperBound.Fill( 0.0 );
	  lowerBound.Fill( 0.0 );

	  optimizer->SetBoundSelection( boundSelect );
	  optimizer->SetUpperBound( upperBound );
	  optimizer->SetLowerBound( lowerBound );
	  optimizer->SetCostFunctionConvergenceFactor( 1e+7 );
	  optimizer->SetGradientConvergenceTolerance( 1.0e-35 );
	  optimizer->SetNumberOfIterations( 30 );
	  optimizer->SetMaximumNumberOfFunctionEvaluations( 200 );
	  optimizer->SetMaximumNumberOfCorrections( 6 );

	  if (imageIndex==1) {
	  	  optimizer->SetInitialPosition(initialTransform1->GetParameters());
	  } else {
	   	  optimizer->SetInitialPosition(transformParameters);
	  }

	  observer = CommandIterationUpdate::New();
	  optimizer->AddObserver( itk::IterationEvent(), observer );

	  command = CommandType::New();
	  registration->AddObserver( itk::MultiResolutionIterationEvent(), command );

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

	  double metricValue1 = optimizer->GetCurrentMetricValue();
	  outfile_metric << metricValue1 << std::endl;

	  // Compare the three metric values after optimization.
	  // Choose the one with the lowest value.
	  //
	  if (imageIndex==1) {
		  if (metricValue1<=metricValue2&&metricValue1<=metricValue3) {
			  outputTransform = initialTransform1->Clone();
		  }
		  if (metricValue2<=metricValue1&&metricValue2<=metricValue3) {
			  outputTransform = initialTransform2->Clone();
		  }
		  if (metricValue3<=metricValue1&&metricValue3<=metricValue2) {
			  outputTransform = initialTransform3->Clone();
		  }
	  } else {
		  if (metricValue1<=metricValue2&&metricValue1<=metricValue3) {
			  outputTransform = outputTransform1->Clone();
		  }
		  if (metricValue2<=metricValue1&&metricValue2<=metricValue3) {
			  outputTransform = outputTransform2->Clone();
		  }
		  if (metricValue3<=metricValue1&&metricValue3<=metricValue2) {
			  outputTransform = outputTransform3->Clone();
		  }
	  }

	  transformParameters = outputTransform->GetParameters();
	  for (int i=0; i<6; i++)
	  {
		  t[i] = transformParameters[i];
	  }
	  metric->SetPreviousTransformParameters(t, 6);
	  std::cout << "Center: " << outputTransform->GetCenter() << std::endl;
	  std::cout << "Matrix: " << outputTransform->GetMatrix() << std::endl;
	  std::cout << "Offset: " << outputTransform->GetOffset() << std::endl;
	  outputTransform->GetInverse(inverseOutputTransform);

      std::cout << "output transform:" << outputTransform->GetParameters() << std::endl;

      outfile_parameters << movedImageName << std::endl;
      outfile_parameters << transformParameters << std::endl;
      outfile_parameters << inverseOutputTransform->GetParameters() << std::endl;

	  // Finally we use the last transform in order to resample the image.
	  //
	  typedef itk::ResampleImageFilter<
								MovingImageType,
								FixedImageType >    ResampleFilterType;

	  ResampleFilterType::Pointer resample = ResampleFilterType::New();

	  resample->SetTransform( outputTransform );
	  resample->SetInput( movingImage );

	  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
	  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
	  resample->SetOutputSpacing( fixedImage->GetSpacing() );
	  resample->SetOutputDirection( fixedImage->GetDirection() );
	  resample->SetDefaultPixelValue( 0 );

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
	  //
	  typedef itk::Vector< float, ImageDimension >          VectorPixelType;
	  typedef itk::Image< VectorPixelType, ImageDimension > DisplacementFieldImageType;

	  typedef itk::TransformToDisplacementFieldFilter<
							DisplacementFieldImageType,
							CoordinateRepType >             DisplacementFieldGeneratorType;

	  DisplacementFieldGeneratorType::Pointer dispfieldGenerator =
													 DisplacementFieldGeneratorType::New();
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

	  fieldWriter->SetInput( dispfieldGenerator->GetOutput() );
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

	  // Finally we use the last inverse transform in order to resample the image.
	  //
	  ResampleFilterType::Pointer resampleInverse = ResampleFilterType::New();

	  resampleInverse->SetTransform( inverseOutputTransform );
	  resampleInverse->SetInput( fixedImage );

	  resampleInverse->SetSize(    movingImage->GetLargestPossibleRegion().GetSize() );
	  resampleInverse->SetOutputOrigin(  movingImage->GetOrigin() );
	  resampleInverse->SetOutputSpacing( movingImage->GetSpacing() );
	  resampleInverse->SetOutputDirection( movingImage->GetDirection() );
	  resampleInverse->SetDefaultPixelValue( 0 );

	  WriterType::Pointer      writerInverse =  WriterType::New();
	  CastFilterType::Pointer  casterInverse =  CastFilterType::New();


	  writerInverse->SetFileName( inverseMovedImageName );


	  casterInverse->SetInput( resampleInverse->GetOutput() );
	  writerInverse->SetInput( casterInverse->GetOutput()   );

	  try
		{
		  writerInverse->Update();
		}
	  catch( itk::ExceptionObject & err )
		{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
		}

	  // Generate the explicit inverse deformation field resulting from
	  // the registration.
	  //
	  DisplacementFieldGeneratorType::Pointer inverseDispfieldGenerator =
													 DisplacementFieldGeneratorType::New();
	  inverseDispfieldGenerator->UseReferenceImageOn();
	  inverseDispfieldGenerator->SetReferenceImage( movingImage );
	  inverseDispfieldGenerator->SetTransform( inverseOutputTransform );
	  try
		{
		  inverseDispfieldGenerator->Update();
		}
	  catch ( itk::ExceptionObject & err )
		{
		std::cerr << "Exception detected while generating deformation field";
		std::cerr << " : "  << err << std::endl;
		return EXIT_FAILURE;
		}

	  FieldWriterType::Pointer inverseFieldWriter = FieldWriterType::New();

	  inverseFieldWriter->SetInput( inverseDispfieldGenerator->GetOutput() );
	  inverseFieldWriter->SetFileName( inverseWarpFieldName );
	  try
	    {
		  inverseFieldWriter->Update();
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
