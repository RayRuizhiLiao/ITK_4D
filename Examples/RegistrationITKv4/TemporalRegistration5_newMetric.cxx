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
#include "itkLBFGSBOptimizerv4.h"
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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"

#include "itkIdentityTransform.h"

#include "itkBSplineTransformInitializer.h"
#include "itkTransformToDisplacementFieldFilter.h"

#include "itkANTSNeighborhoodCorrelationImageToImageTemporalMetricv4.h"


//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  typedef itk::LBFGSBOptimizerv4     OptimizerType;
  typedef   const OptimizerType *    OptimizerPointer;

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
  std::cout << optimizer->GetCurrentMetricValue() << "   ";
  std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
  }
};

template <typename T>
std::string to_string(T value)
{
	std::ostringstream os;
	os << value;
	return os.str();
}


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
  const unsigned int SpaceDimension = ImageDimension;
  const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;

  typedef itk::BSplineTransform<
                            CoordinateRepType,
                            SpaceDimension,
                            SplineOrder >     TransformType;
  // Software Guide : EndCodeSnippet


  typedef itk::LBFGSBOptimizerv4       OptimizerType;


  //typedef itk::MeanSquaresImageToImageMetricv4<
  //                                  FixedImageType,
  //                                  MovingImageType >    MetricType;

  typedef itk::ANTSNeighborhoodCorrelationImageToImageTemporalMetricv4<
                                          FixedImageType,
                                          MovingImageType >    MetricType;

  typedef itk::ImageRegistrationMethodv4<
                                    FixedImageType,
                                    MovingImageType >    RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();


  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );

  int numOfImages  = 0;
  std::sscanf(argv[2], "%d", &numOfImages);
  std::cout << "There are " << numOfImages << " images in the given series!" << std::endl;


  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  OptimizerType::ParametersType transformParameters;

  std::string outputFolder(argv[1]);

  for (int imageIndex=1; imageIndex<numOfImages; imageIndex++) {

	  fixedImageReader->SetFileName(  argv[3] );
	  movingImageReader->SetFileName( argv[3+imageIndex] );

	  std::string fixedImageName(argv[3]);
	  std::string slash = "/";
	  std::size_t fixedSlashIndex = fixedImageName.find(slash, fixedImageName.length()-25);
	  std::string movingImageName(argv[3+imageIndex]);
	  std::size_t movingSlashIndex = movingImageName.find(slash, movingImageName.length()-25);

	  std::string movedImageName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + ".nii.gz";
	  std::cout << "Moved Image Name: " << movedImageName << std::endl;
	  std::string warpFieldName = outputFolder + movingImageName.substr(movingSlashIndex+1, movingImageName.length()-movingSlashIndex-8) + "_to_" + fixedImageName.substr(fixedSlashIndex+1, fixedImageName.length()-fixedSlashIndex-8) + "_warp.nii.gz";
	  std::cout << "Warp Image Name: " << warpFieldName << std::endl;

	  FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();

	  registration->SetFixedImage(  fixedImage   );
	  registration->SetMovingImage(   movingImageReader->GetOutput()   );

	  fixedImageReader->Update();

	  //  Software Guide : BeginLatex
	  //
	  //  The transform object is constructed, initialized like previous examples
	  //  and passed to the registration method.
	  //
	  //  \index{itk::ImageRegistrationMethodv4!SetInitialTransform()}
	  //  \index{itk::ImageRegistrationMethodv4!InPlaceOn()}
	  //
	  //  Software Guide : EndLatex

	  // Software Guide : BeginCodeSnippet
	  TransformType::Pointer  outputBSplineTransform = TransformType::New();
	  // Software Guide : EndCodeSnippet

	  // Initialize the transform
	  typedef itk::BSplineTransformInitializer< TransformType,
												FixedImageType>      InitializerType;

	  InitializerType::Pointer transformInitializer = InitializerType::New();

	  //unsigned int numberOfGridNodesInOneDimension = 8;

	  TransformType::MeshSizeType             meshSize;
	  meshSize.SetElement(0, 20);
	  meshSize.SetElement(1, 20);
	  meshSize.SetElement(2, 16);
	  //meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );

	  transformInitializer->SetTransform( outputBSplineTransform );
	  transformInitializer->SetImage( fixedImage );
	  transformInitializer->SetTransformDomainMeshSize( meshSize );
	  transformInitializer->InitializeTransform();

	  // Set transform to identity
	  typedef TransformType::ParametersType     ParametersType;
	  const unsigned int numberOfParameters =
				   outputBSplineTransform->GetNumberOfParameters();
	  ParametersType parameters( numberOfParameters );
	  parameters.Fill( 0.0 );
	  outputBSplineTransform->SetParameters( parameters );
	  TransformType::Pointer  identityBSplineTransform = outputBSplineTransform->Clone();

	  // Software Guide : BeginCodeSnippet
	  registration->SetInitialTransform( outputBSplineTransform );
	  registration->InPlaceOn();
	  // Software Guide : EndCodeSnippet

	  //  A single level registration process is run using
	  //  the shrink factor 1 and smoothing sigma 0.
	  //
	  const unsigned int numberOfLevels = 1;

	  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	  shrinkFactorsPerLevel.SetSize( numberOfLevels );
	  shrinkFactorsPerLevel[0] = 1;
	  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	  smoothingSigmasPerLevel.SetSize( numberOfLevels );
	  smoothingSigmasPerLevel[0] = 0;

	  registration->SetNumberOfLevels( numberOfLevels );
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
		outputBSplineTransform->GetNumberOfParameters();
	  OptimizerType::BoundSelectionType boundSelect( numParameters );
	  OptimizerType::BoundValueType upperBound( numParameters );
	  OptimizerType::BoundValueType lowerBound( numParameters );

	  boundSelect.Fill( OptimizerType::UNBOUNDED );
	  upperBound.Fill( 0.0 );
	  lowerBound.Fill( 0.0 );

	  optimizer->SetBoundSelection( boundSelect );
	  optimizer->SetUpperBound( upperBound );
	  optimizer->SetLowerBound( lowerBound );
	  optimizer->TraceOn();
	  optimizer->SetCostFunctionConvergenceFactor( 1e+11 );
	  optimizer->SetGradientConvergenceTolerance( 1.0e-10 );
	  optimizer->SetNumberOfIterations( 30 );
	  optimizer->SetMaximumNumberOfFunctionEvaluations( 50 );
	  optimizer->SetMaximumNumberOfCorrections( 5 );
	  std::cout << "Number of Threads: " << optimizer->GetNumberOfThreads() << std::endl;
	  if (imageIndex==1) {
		  optimizer->SetInitialPosition(outputBSplineTransform->GetParameters());
	  } else {
		  optimizer->SetInitialPosition(transformParameters);
	  }
	  //std::cout << "Initialization: " << optimizer->GetInitialPosition() << std::endl;
	  // Software Guide : EndCodeSnippet

	  // Create the Command observer and register it with the optimizer.
	  //
	  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	  optimizer->AddObserver( itk::IterationEvent(), observer );


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
	  transformParameters = outputBSplineTransform->GetParameters();
	  // Software Guide : EndCodeSnippet

	  std::cout << "Last Transform Parameters" << std::endl;
	  //std::cout << transformParameters << std::endl;

	  // Finally we use the last transform in order to resample the image.
	  //
	  typedef itk::ResampleImageFilter<
								MovingImageType,
								FixedImageType >    ResampleFilterType;

	  ResampleFilterType::Pointer resample = ResampleFilterType::New();

	  resample->SetTransform( outputBSplineTransform );
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
	  dispfieldGenerator->UseReferenceImageOn();
	  dispfieldGenerator->SetReferenceImage( fixedImage );
	  dispfieldGenerator->SetTransform( outputBSplineTransform );
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

  }

  return EXIT_SUCCESS;
}
