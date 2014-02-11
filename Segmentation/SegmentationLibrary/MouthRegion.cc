//MouthRegion.cc
//Author - Bud Goswami
//Defines the functions that use OpenCV and Omni to detect the mouth region

#include "MouthRegion.hh"
// Create memory for calculations

Point2dC topleft(0,0);
Point2dC botright(0,0);
// Function prototype for detecting and drawing an object from an image
void detect_and_draw( IplImage* img, CvHaarClassifierCascade* cascade )
{
	int scale = 1;
	CvMemStorage* storage = cvCreateMemStorage(0);
	 // Create a new image based on the input image
	 IplImage* temp = cvCreateImage( cvSize(img->width/scale,img->height/scale), 8, 3 );

	 // Create two points to represent the face locations
	 CvPoint pt1, pt2;
	 int i;
	 // Clear the memory storage which was used before
	 cvClearMemStorage(storage);
	 // Find whether the cascade is loaded, to find the faces. If yes, then:
	 if( cascade )
	 {
		 	// There can be more than one face in an image. So create a growable sequence of faces.
	        // Detect the objects and store them in the sequence
	        CvSeq* faces = cvHaarDetectObjects( img, cascade, storage,1.1, 2, CV_HAAR_DO_CANNY_PRUNING,cvSize(40, 40) );
	        // Loop the number of faces found.
	        for( i = 0; i < (faces ? faces->total : 0); i++ )
	        {
	        	// Create a new rectangle for drawing the face
	            CvRect* r = (CvRect*)cvGetSeqElem( faces, i );
	            // Find the dimensions of the face,and scale it if necessary
	            pt1.x = r->x*scale;
	            pt2.x = (r->x+r->width)*scale;
	            pt1.y = r->y*scale;
	            pt2.y = (r->y+r->height)*scale;
	            topleft.Row() = pt1.y;
	            topleft.Col() = pt1.x;
	            botright.Row() = pt2.y;
	            botright.Col() = pt2.x;            
	            // Draw the rectangle in the input image
	            cvRectangle( img, pt1, pt2, CV_RGB(255,0,0), 3, 8, 0 );
	        }
	    }
	    // Release the temp image created.
	    cvReleaseImage( &temp );
}

//Function definition to obtain Face Co-ordinates
ImageRectangleC GetFaceCoords(const FilenameC &img_name, const FilenameC &c_name)
{
		IplImage* image = 0;
		image = cvLoadImage( img_name.chars(), 1 );
		//IMPORTANT - GOING THROUGH RAVLIPLIMAGE CAUSES ISSUES AND OPENCV SAYS UNSUPPORTED IMAGE TYPE
		//Go through RAVL
		//ImageC<RealRGBValueC> rgb_img;
		//if(!Load(img_name, rgb_img)) cerr<<"Could not load RGB Image"<<endl;
		//if(!RavlImage2IplImage(rgb_img,image) ) cerr<<"Image conversion to OpenCV not possible"<<endl;
		static CvHaarClassifierCascade* cascade = 0;
		cascade = (CvHaarClassifierCascade*)cvLoad( c_name.chars(), 0, 0, 0 ); 
		// Check whether the cascade has loaded successfully. Else report and error and quit
		if( !cascade )
		{
			fprintf( stderr, "ERROR: Could not load classifier cascade\n" );
			//return ImageRectangleC(0,0);
		} 
		//static CvMemStorage* storage = 0;
		//storage = cvCreateMemStorage(0);
		if( image )
		{
			//Detect and draw the face
			detect_and_draw( image, cascade );
			// Release the image memory
			cvReleaseImage( &image );	
   	 }
 		ImageRectangleC facerec(topleft, botright);
		return facerec;
}
ImageRectangleC GetFaceCoords(const ImageC<RealRGBValueC> &img, const FilenameC &c_name)
{
		IplImage* image = 0;
		RavlImage2IplImage(img,image);  
		static CvHaarClassifierCascade* cascade = 0;
		cascade = (CvHaarClassifierCascade*)cvLoad( c_name.chars(), 0, 0, 0 ); 
		// Check whether the cascade has loaded successfully. Else report and error and quit
		if( !cascade )
		{
			fprintf( stderr, "ERROR: Could not load classifier cascade\n" );
			//return ImageRectangleC(0,0);
		} 
		//static CvMemStorage* storage = 0;
		//storage = cvCreateMemStorage(0);
		if( image )
		{
			//Detect and draw the face
			detect_and_draw( image, cascade );
			// Release the image memory
			cvReleaseImage( &image );	
   	 }
 		ImageRectangleC facerec(topleft, botright);
		return facerec;
}
Tuple2C<Point2dC,Point2dC> GetEyeLoc(const ImageC<RealRGBValueC> &img)
{
	//FilenameC  faceDetectionModelHQ = "/home/bud/WorkHome/lip-tracking/OmniFaceModels/HighQualityFD.abs";
	//FilenameC  faceDetectionModelLQ = "/home/bud/WorkHome/lip-tracking/OmniFaceModels/LowQualityFaceScan.abs";
	StringC OMNIHOME = "/vol/vssp/lip-tracking/OMNIModels/";
	FilenameC  faceDetectionModelHQ =FilenameC(OMNIHOME + "HighQualityFD.abs");
	FilenameC  faceDetectionModelLQ =FilenameC(OMNIHOME + "LowQualityFaceScan.abs");
	//FilenameC  faceDetectionModelHQ =OMNIHOME"/share/FaceDetect/models/HighQualityFD.abs";
	//FilenameC  faceDetectionModelLQ =OMNIHOME"/share/FaceDetect/models/LowQualityFaceScan.abs";
	//FilenameC  faceDetectionModelHQ = "/vol/vssp/localsoft/External/FaceDetect/models/HighQualityFD.abs";
	//FilenameC  faceDetectionModelLQ = "/vol/vssp/localsoft/External/FaceDetect/models/LowQualityFaceScan.abs";
	FaceScanC detectAAMFace;           // face locator that can provide face shape from AAM
	DetectFaceC detectFace;            // eezy-use face locator
	DetectedFaceC detectedFace;
	//cout<<"set up face and eye locators"<<endl;
   FaceDetectionModelC faceDetectionModel;
   if(!Load(faceDetectionModelHQ, faceDetectionModel))
   	RavlIssueError("Trouble loading eezy-use face detection model\n");
   detectFace = DetectFaceC(faceDetectionModel,1.); 
   if(!Load(faceDetectionModelLQ, detectAAMFace))
    	RavlIssueError("Trouble loading AAM face detection model\n");
   ImageC<ByteRGBValueC> fullSizeImg = RealRGBImageCT2ByteRGBImageCT(img);
   UIntT factor(1);
   //cout<<"DetectFace.Apply()"<<endl;
   DetectedFaceC tmp = detectFace.Apply(fullSizeImg);
   tmp.LeftEye()  *= factor;
   tmp.RightEye() *= factor;
   //cout<<"About to Get Co-ords"<<endl;
   detectedFace = DetectedFaceC(fullSizeImg, tmp.LeftEye(), tmp.RightEye(), tmp.Quality());
   //cout<<"Convert to output format"<<endl;
   Tuple2C<Point2dC,Point2dC> eyes(Point2dC(detectedFace.LeftEye()[0], detectedFace.LeftEye()[1]),Point2dC(detectedFace.RightEye()[0],detectedFace.RightEye()[1])); 		
   return eyes;
}

Tuple2C<Point2dC,Point2dC> GetEyeLoc(const FilenameC &ImageFile)
{
	ImageC<RealRGBValueC> fullSizeImg;
    RavlAlwaysAssertMsg(Load(ImageFile, fullSizeImg), "Could not load image: " + ImageFile);	
    Tuple2C<Point2dC,Point2dC> eyes = GetEyeLoc(fullSizeImg);
    return eyes;
}
ImageRectangleC GetMouthRegion(const FilenameC &img_name, const FilenameC &c_name, const IndexC &window = (IndexC)10)
{
	ImageRectangleC facerec = GetFaceCoords(img_name,c_name);
	ImageRectangleC res;
	try
	{
		Tuple2C<Point2dC,Point2dC> eyes = GetEyeLoc(img_name);
		IndexC trow = facerec.BRow() - (facerec.Rows())/3;
		IndexC mincol = eyes.Data2().Col() - ((eyes.Data2().Col() - facerec.LCol())/4);
		IndexC maxcol = eyes.Data1().Col() + ((facerec.RCol() - eyes.Data1().Col())/4);
		ImageRectangleC botfacehalf(trow,facerec.BRow(),mincol, maxcol);
		res = botfacehalf;
	}
	catch (OmniN::ColossusExceptionNoFaceFoundC er)
  	{
   	//If eye estimation doesn't work, just return the bottom half of the face clipped by a window
   	ImageRectangleC bothalf((facerec.TRow() + facerec.BRow())/2,facerec.BRow() - window,facerec.LCol() + (2*window), facerec.RCol() - (2*window));	 	
   	res = bothalf;
  	}
	return res;
}
ImageRectangleC GetMouthRegion(const ImageC<RealRGBValueC> &img, const FilenameC &c_name, const IndexC &window = (IndexC)10)
{
	ImageRectangleC facerec = GetFaceCoords(img,c_name);
	ImageRectangleC res;
	try
	{
		Tuple2C<Point2dC,Point2dC> eyes = GetEyeLoc(img);
		IndexC trow = facerec.BRow() - (facerec.Rows())/3;
		IndexC mincol = eyes.Data2().Col() - ((eyes.Data2().Col() - facerec.LCol())/4);
		IndexC maxcol = eyes.Data1().Col() + ((facerec.RCol() - eyes.Data1().Col())/4);
		ImageRectangleC botfacehalf(trow,facerec.BRow(),mincol, maxcol);
		res = botfacehalf;
	}
	catch (OmniN::ColossusExceptionNoFaceFoundC er)
  	{
   	//If eye estimation doesn't work, just return the bottom half of the face clipped by a window
   	ImageRectangleC bothalf((facerec.TRow() + facerec.BRow())/2,facerec.BRow() - window,facerec.LCol() + (2*window), facerec.RCol() - (2*window));	 	
   	res = bothalf;
  	}
	return res;
}
