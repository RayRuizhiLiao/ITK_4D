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

#include "itkRegionalMaximaImageFilter.h"
#include "itkHConvexImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkTestingMacros.h"


template< typename TInputImage, typename TOutputImage >
int RegionalMaximaImageFilterTestHelper( std::string inputImageFile,
  std::string outputImageFile, std::string outputImageFile2,
  bool fullyConnected, bool flatIsMaxima )
{
  typedef TInputImage InputImageType;
  typedef TInputImage OutputImageType;

  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFile );

  typedef itk::RegionalMaximaImageFilter< InputImageType, OutputImageType >
    FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( reader->GetOutput() );

  TEST_SET_GET_BOOLEAN( filter, FullyConnected, fullyConnected );

  TEST_SET_GET_BOOLEAN( filter, FlatIsMaxima, flatIsMaxima );

  typename FilterType::OutputImagePixelType foregroundValue =
    itk::NumericTraits< typename FilterType::OutputImagePixelType >::max();
  filter->SetForegroundValue( foregroundValue );
  TEST_SET_GET_VALUE( foregroundValue, filter->GetForegroundValue() );

  typename FilterType::OutputImagePixelType backgroundValue =
    itk::NumericTraits< typename FilterType::OutputImagePixelType >::NonpositiveMin();
  filter->SetBackgroundValue( backgroundValue );
  TEST_SET_GET_VALUE( backgroundValue, filter->GetBackgroundValue() );

  itk::SimpleFilterWatcher watcher( filter, "RegionalMaximaImageFilter" );

  // Write the output images
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( outputImageFile );
  writer->Update();


  // Produce the same output with other filters
  typedef itk::HConvexImageFilter< InputImageType, InputImageType >
    ConvexFilterType;
  typename ConvexFilterType::Pointer convexFilter = ConvexFilterType::New();
  convexFilter->SetInput( reader->GetOutput() );
  convexFilter->SetFullyConnected( fullyConnected );
  convexFilter->SetHeight( 1 );

  // Convex gives maxima with value = 1 and others with value = 0
  // Rescale the image so we have maxima = 255 other = 0
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType >
    RescaleFilterType;
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetInput( convexFilter->GetOutput() );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetOutputMinimum( 0 );

  typename WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( rescaler->GetOutput() );
  writer2->SetFileName( outputImageFile2 );
  writer2->Update();


  std::cerr << "Test finished" << std::endl;
  return EXIT_SUCCESS;
}

// A test routine for regional extrema using flooding
int itkRegionalMaximaImageFilterTest( int argc, char * argv[] )
{
  if( argc < 7 )
    {
    std::cerr << "Missing parameters" << std::endl;
    std::cerr << "Usage: " << argv[0]
      << " inputImageFile"
      << " outputImageFile"
      << " outputImageFile2"
      << " dimension"
      << " fullyConnected"
      << " flatIsMaxima";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputImageFile = argv[1];
  std::string outputImageFile = argv[2];
  std::string outputImageFile2 = argv[3];

  unsigned int dimension = atoi( argv[4] );

  bool fullyConnected = atoi( argv[5] );
  bool flatIsMaxima = atoi( argv[6] );

  // Exercise basic object methods
  // Done outside the helper function in the test because GCC is limited
  // when calling overloaded base class functions.
  typedef unsigned char               PixelType;
  typedef itk::Image< PixelType, 2 >  ImageType;

  typedef itk::RegionalMaximaImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();

  EXERCISE_BASIC_OBJECT_METHODS( filter, RegionalMaximaImageFilter,
    ImageToImageFilter );

  if( dimension == 2 )
    {
    typedef itk::Image< PixelType, 2 > Image2DType;
    return RegionalMaximaImageFilterTestHelper< Image2DType, Image2DType >(
      inputImageFile, outputImageFile, outputImageFile2, fullyConnected,
      flatIsMaxima );
    }
  else if( dimension == 3 )
    {
    typedef itk::Image< PixelType, 3 > Image3DType;
    return RegionalMaximaImageFilterTestHelper< Image3DType, Image3DType >(
      inputImageFile, outputImageFile, outputImageFile2, fullyConnected,
      flatIsMaxima );
    }
  else
    {
    std::cerr << "Test failed!" << std::endl;
    std::cerr << "Unsupported dimension: " << dimension << std::endl;
    std::cerr << "Only dimensions 2 and 3 are supported." << std::endl;
    return EXIT_FAILURE;
    }
}
