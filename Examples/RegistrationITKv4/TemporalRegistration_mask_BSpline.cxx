/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

// Software Guide : BeginLatex
//
// This example illustrates the use of the \doxygen{BSplineTransform}
// class for performing registration of two $3D$ images. The example code is
// for the most part identical to the code presented in
// Section~\ref{sec:BSplinesMultiGridImageRegistration}. The major difference is
// that in this example we set the image dimension to 3 and replace the
// \doxygen{LBFGSOptimizerv4} optimizer with the \doxygen{LBFGSBOptimizerv4}. We
// made the modification because we found that LBFGS does not behave well when
// the starting position is at or close to optimal; instead we used LBFGSB in
// unconstrained mode.
//
//
// \index{itk::BSplineTransform}
// \index{itk::BSplineTransform!DeformableRegistration}
// \index{itk::LBFGSBOptimizerv4}
//
//
// Software Guide : EndLatex

#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

//  Software Guide : BeginLatex
//
//  The following are the most relevant headers to this example.
//
//  \index{itk::BSplineTransform!header}
//  \index{itk::LBFGSBOptimizerv4!header}
//
//  Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
#include "itkBSplineTransform.h"
#include "itkVersorRigid3DTransform.h"

#include "itkLBFGSBOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"

// Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//
//  The parameter space of the \code{BSplineTransform} is composed by
//  the set of all the deformations associated with the nodes of the BSpline
//  grid.  This large number of parameters enables it to represent a wide
//  variety of deformations, at the cost of requiring a
//  significant amount of computation time.
//
//  \index{itk::BSplineTransform!header}
//
//  Software Guide : EndLatex
#include "itkSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"

#include "itkIdentityTransform.h"

#include "itkBSplineTransformInitializer.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkInverseDisplacementFieldImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkIterativeInverseDisplacementFieldImageFilter.h"
#include "itkInverseDeformationFieldImageFilter.h"


#include "itkANTSNeighborhoodCorrelationImageToImageTemporalMetricv4.h"


//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"

// Software Guide : BeginCodeSnippet
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // We then define \code{Self}, \code{Superclass}, \code{Pointer},
  // \code{New()} and a constructor in a similar fashion to the
  // \code{CommandIterationUpdate} class in Section
  // \ref{sec:MonitoringImageRegistration}.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );

protected:
  RegistrationInterfaceCommand() {};
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // For convenience, we declare types useful for converting pointers
  // in the \code{Execute()} method.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
public:
  typedef   TRegistration      RegistrationType;
  typedef   RegistrationType * RegistrationPointer;
  typedef itk::LBFGSBOptimizerv4       OptimizerType;
  typedef   OptimizerType * OptimizerPointer;
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // Two arguments are passed to the \code{Execute()} method: the first
  // is the pointer to the object which invoked the event and the
  // second is the event that was invoked.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  void Execute( itk::Object * object,
                const itk::EventObject & event) ITK_OVERRIDE
    {
    // Software Guide : EndCodeSnippet

    // Software Guide : BeginLatex
    //
    // First we verify that the event invoked is of the right type,
    // \code{itk::MultiResolutionIterationEvent()}.
    // If not, we return without any further action.
    //
    // Software Guide : EndLatex

    // Software Guide : BeginCodeSnippet
    if( !(itk::MultiResolutionIterationEvent().CheckEvent( &event ) ) )
      {
      return;
      }
    // Software Guide : EndCodeSnippet

    // Software Guide : BeginLatex
    //
    // We then convert the input object pointer to a RegistrationPointer.
    // Note that no error checking is done here to verify the
    // \code{dynamic\_cast} was successful since we know the actual object
    // is a registration method. Then we ask for the optimizer object
    // from the registration method.
    //
    // Software Guide : EndLatex

    // Software Guide : BeginCodeSnippet
    RegistrationPointer registration =
      static_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer =  static_cast< OptimizerPointer >(
        registration->GetModifiableOptimizer() );
    // Software Guide : EndCodeSnippet

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
    std::cout << optimizer->GetValue() << std::endl;
    //std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
    catch(itk::ExceptionObject)
    {
    }

    // Software Guide : BeginLatex
    //
    // If this is the first resolution level we set the learning rate
    // (representing the first step size) and the minimum step length (representing
    // the convergence criterion) to large values.  At each subsequent resolution
    // level, we will reduce the minimum step length by a factor of 5 in order to
    // allow the optimizer to focus on progressively smaller regions. The learning
    // rate is set up to the current step length. In this way, when the
    // optimizer is reinitialized at the beginning of the registration process for
    // the next level, the step length will simply start with the last value used
    // for the previous level. This will guarantee the continuity of the path
    // taken by the optimizer through the parameter space.
    //
    // Software Guide : EndLatex

    // Software Guide : BeginCodeSnippet
    if ( registration->GetCurrentLevel() == 0 )
      {
      //optimizer->SetLearningRate( 5.00 );
      //optimizer->SetLearningRate( 1 );
      //optimizer->SetMinimumStepLength( 0.5 );
      }
    else
      {
    	optimizer->SetInitialPosition(optimizer->GetCurrentPosition());
      //optimizer->SetLearningRate( optimizer->GetCurrentStepLength() );
      //optimizer->SetMinimumStepLength(
      //  optimizer->GetMinimumStepLength() * 0.2 );
      }
    // Software Guide : EndCodeSnippet
    }

  // Software Guide : BeginLatex
  //
  // Another version of the \code{Execute()} method accepting a \code{const}
  // input object is also required since this method is defined as pure virtual
  // in the base class.  This version simply returns without taking any action.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
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
  void SetOutputOptimizationLog (std::string outputFile) {
	  this->outputOptimizationLog = outputFile;
  }
  std::string outputOptimizationLog;

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
  {
  Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
  {
  std::ofstream outfile2;
  outfile2.open(this->outputOptimizationLog.c_str(), std::ofstream::out|std::ofstream::app);

  OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
  if( !(itk::IterationEvent().CheckEvent( &event )) )
    {
    return;
    }
  std::cout << optimizer->GetCurrentIteration() << "   ";
  std::cout << optimizer->GetValue() << "   ";
 // std::cout << optimizer->GetCurrentPosition() << "   ";
  //std::cout << optimizer->GetCachedDerivative() << "   ";
  std::cout << m_CumulativeIterationIndex++ << std::endl;
  outfile2.close();
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


  //  Software Guide : BeginLatex
  //
  //  We instantiate now the type of the \code{BSplineTransform} using
  //  as template parameters the type for coordinates representation, the
  //  dimension of the space, and the order of the BSpline.
  //
  //  \index{BSplineTransform!New}
  //  \index{BSplineTransform!Instantiation}
  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef double CoordinateRepType;

  const unsigned int SpaceDimension = ImageDimension;
  const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;

  typedef itk::BSplineTransform<
                            CoordinateRepType,
                            SpaceDimension,
                            SplineOrder >     TransformType;
  // Software Guide : EndCodeSnippet


  typedef itk::LBFGSBOptimizerv4       OptimizerType;

  typedef itk::ANTSNeighborhoodCorrelationImageToImageTemporalMetricv4<
                                          FixedImageType,
                                          MovingImageType >    MetricType;

  typedef itk::ImageRegistrationMethodv4<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

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

  std::string outputFolder(argv[1]);

  maskImageReader->SetFileName( argv[3] );

  std::ofstream outfile;
  std::string outputMetricValues = outputFolder + "metric_values.txt";
  outfile.open(outputMetricValues.c_str(), std::ofstream::out|std::ofstream::app);

  std::ofstream outfile3;
  std::string outputTransformParameters = outputFolder + "transform_parameters.txt";
  outfile3.open(outputTransformParameters.c_str(), std::ofstream::out|std::ofstream::app);

  std::ofstream outfile4;
  std::string outputValidPoints = outputFolder + "valid_points.txt";
  outfile4.open(outputValidPoints.c_str(), std::ofstream::out|std::ofstream::app);

  double w1;
  std::sscanf(argv[4], "%lf", &w1);
  double w2;
  std::sscanf(argv[5], "%lf", &w2);

  double* t = new double[6];
  for (int i=0; i<6; i++)
  {
	  t[i] = 0;
  }

  // Software Guide : BeginCodeSnippet
  TransformType::Pointer  outputTransform = TransformType::New();
  TransformType::Pointer  inverseOutputTransform = TransformType::New();
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
	  std::string warpFieldName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_warp.nii.gz";
	  std::cout << "Warp Image Name: " << warpFieldName << std::endl;

	  std::string inverseMovedImageName = outputFolder + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_to_" + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + ".nii.gz";
	  std::cout << "Inverse Moved Image Name: " << inverseMovedImageName << std::endl;
	  std::string inverseWarpFieldName = outputFolder + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_to_" + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_inverseWarp.nii.gz";
	  std::cout << "Inverse Warp Image Name: " << inverseWarpFieldName << std::endl;

	  std::ofstream outfile2;
	  std::string outputOptimizationLog = outputFolder + "optimization_log.txt";
	  outfile2.open(outputOptimizationLog.c_str(), std::ofstream::out|std::ofstream::app);
      outfile2 << movedImageName << std::endl;
      outfile2.close();

      MetricType::Pointer         metric        = MetricType::New();
      OptimizerType::Pointer      optimizer     = OptimizerType::New();
      RegistrationType::Pointer   registration  = RegistrationType::New();

      registration->SetMetric(        metric        );
      registration->SetOptimizer(     optimizer     );

      metric->SetTemporalSmoothness1(w1);
      metric->SetTemporalSmoothness2(w2);
      metric->SetPreviousTransformParameters(t, 6);

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


	  //  Software Guide : BeginLatex
	  //
	  //  The transform object is constructed, initialized like previous examples
	  //  and passed to the registration method.
	  //
	  //  \index{itk::ImageRegistrationMethodv4!SetInitialTransform()}
	  //  \index{itk::ImageRegistrationMethodv4!InPlaceOn()}
	  //
	  //  Software Guide : EndLatex


	  outfile3 << outputTransform->GetParameters() << std::endl;
	  outfile4 << metric->GetPreviousNumberOfValidPoints() << std::endl;

	  // Initialize the transform
	  typedef itk::BSplineTransformInitializer< TransformType,
												FixedImageType>      InitializerType;
	  InitializerType::Pointer transformInitializer = InitializerType::New();

	  //unsigned int numberOfGridNodesInOneDimension = 8;

	  TransformType::MeshSizeType             meshSize;
	  meshSize.SetElement(0, 7);
	  meshSize.SetElement(1, 7);
	  meshSize.SetElement(2, 5);
	  //meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );

	  if (imageIndex==1) {
		  transformInitializer->SetTransform(   initialTransform );
		  transformInitializer->SetImage(  fixedImageReader->GetOutput() );
		  transformInitializer->SetTransformDomainMeshSize( meshSize );
		  transformInitializer->InitializeTransform();

		  registration->SetInitialTransform( initialTransform );
		  registration->InPlaceOn();
	  } else {
		  transformInitializer->SetTransform(   outputTransform );
		  transformInitializer->SetImage(  fixedImageReader->GetOutput() );
		  transformInitializer->SetTransformDomainMeshSize( meshSize );
		  transformInitializer->InitializeTransform();

		  outfile3 << outputTransform->GetParameters() << std::endl;

		  registration->SetInitialTransform( outputTransform );
	  }

	  //  A single level registration process is run using
	  //  the shrink factor 1 and smoothing sigma 0.
	  //
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

	  //  Software Guide : BeginLatex
	  //
	  //  Next we set the parameters of the LBFGSB Optimizer. Note that
	  //  this optimizer does not support scales estimator and sets all
	  //  the parameters scales to one.
	  //  Also, we should set the boundary condition for each variable, where
	  //  \code{boundSelect[i]} can be set as: \code{UNBOUNDED},
	  //  \code{LOWERBOUNDED}, \code{BOTHBOUNDED}, \code{UPPERBOUNDED}.
	  //
	  //  Software Guide : EndLatex

	  // Software Guide : BeginCodeSnippet
	  const unsigned int numParameters =
			  initialTransform->GetNumberOfParameters();

	  std::cout << "num of transform parameters:" << numParameters << std::endl;


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
	  optimizer->SetMaximumNumberOfCorrections( 5 );
	  //optimizer->SetNumberOfThreads(1);
	  //registration->SetNumberOfThreads(1);
	  //metric->SetMaximumNumberOfThreads(1);
	  if (imageIndex==1) {
	  	  optimizer->SetInitialPosition(initialTransform->GetParameters());
	  } else {
	   	  optimizer->SetInitialPosition(outputTransform->GetParameters());
	  }

	  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	  observer->SetOutputOptimizationLog(outputOptimizationLog);
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
	  outfile << metricValue <<std::endl;

	  // Report the time and memory taken by the registration
	  //chronometer.Report( std::cout );
	  //memorymeter.Report( std::cout );

	  //  Software Guide : BeginLatex
	  //
	  //  Let's execute this example using the rat lung images from the previous examples.
	  //
	  // \begin{itemize}
	  // \item \code{RatLungSlice1.mha}
	  // \item \code{RatLungSlice2.mha}
	  // \end{itemize}
	  //
	  //  The \emph{transform} object is updated during the registration process
	  //  and is passed to the resampler to map the moving image space onto the
	  //  fixed image space.
	  //
	  //  Software Guide : EndLatex

	  // Software Guide : BeginCodeSnippet

	  transformParameters = outputTransform->GetParameters();

      outfile3 << movedImageName << std::endl;
      outfile3 << transformParameters << std::endl;
      outfile3 << inverseOutputTransform->GetParameters() << std::endl;

	  for (int i=0; i<6; i++)
	  {
		  t[i] = transformParameters[i];
	  }

	  // Software Guide : EndCodeSnippet

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

	  DisplacementFieldImageType::Pointer dispField = DisplacementFieldImageType::New();

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

	  dispField = dispfieldGenerator->GetOutput();

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

	  // Finally we use the last inverse transform in order to resample the image.
	  //

	  typedef itk::InverseDeformationFieldImageFilter<
							DisplacementFieldImageType,
							DisplacementFieldImageType >             InverseDisplacementFieldFilterType;

	  InverseDisplacementFieldFilterType::Pointer inverseDispFieldGenerator = InverseDisplacementFieldFilterType::New();

	  inverseDispFieldGenerator->SetInput(dispField);
	  //inverseDispFieldGenerator->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	  //inverseDispFieldGenerator->SetOutputSpacing(fixedImage->GetSpacing());
	  //inverseDispFieldGenerator->SetOutputOrigin(dispfieldGenerator->GetOutputOrigin());
	  inverseDispFieldGenerator->SetSubsamplingFactor(8);
	  inverseDispFieldGenerator->Update();

	  typedef itk::WarpImageFilter<
			  	  	  FixedImageType,
					  MovingImageType,
					  DisplacementFieldImageType >             WarpImageFilterType;

	  WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();

	  warpImageFilter->SetInput(fixedImageReader->GetOutput());
	  warpImageFilter->SetOutputParametersFromImage(fixedImageReader->GetOutput());
	  warpImageFilter->SetDisplacementField(inverseDispFieldGenerator->GetOutput());
	  warpImageFilter->Update();

	  WriterType::Pointer      writer2 =  WriterType::New();

	  writer2->SetFileName( inverseMovedImageName );
	  writer2->SetInput( warpImageFilter->GetOutput() );

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

	  // Generate the explicit inverse deformation field resulting from
	  // the registration.

	  typedef itk::ImageFileWriter< DisplacementFieldImageType >  FieldWriterType;
	  FieldWriterType::Pointer fieldWriter2 = FieldWriterType::New();

	  fieldWriter2->SetInput( inverseDispFieldGenerator->GetOutput() );
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

  outfile.close();
  outfile3.close();
  outfile4.close();

  return EXIT_SUCCESS;
}
