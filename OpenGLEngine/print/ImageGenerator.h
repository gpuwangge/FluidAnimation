#pragma once
#include "app.h"

bool SaveBMP(string filename, int w, int h) {
	// we will store the image data here
	unsigned char *pixels;
	// we get the width/height of the screen into this array

	int screenStats[4];

	// get the width/height of the window
	glGetIntegerv(GL_VIEWPORT, screenStats);

	// generate an array large enough to hold the pixel data 
	// (width*height*bytesPerPixel)
	pixels = new unsigned char[screenStats[2] * screenStats[3] * 3];
	// read in the pixel data, TGA's pixels are BGR aligned
	glReadPixels(0, 0, screenStats[2], screenStats[3], GL_BGR,
		GL_UNSIGNED_BYTE, pixels);

	//Create a new file for writing
	FILE *pFile;
	fopen_s(&pFile, filename.c_str(), "wb");
	if (pFile == NULL) return false;

	BITMAPINFOHEADER BMIH;
	BMIH.biSize = sizeof(BITMAPINFOHEADER);
	BMIH.biSizeImage = w * h * 3;
	// Create the bitmap for this OpenGL context

	BMIH.biSize = sizeof(BITMAPINFOHEADER);
	BMIH.biWidth = w;
	BMIH.biHeight = h;
	BMIH.biPlanes = 1;
	BMIH.biBitCount = 24;
	BMIH.biCompression = BI_RGB;
	BMIH.biSizeImage = w * h * 3;

	BITMAPFILEHEADER bmfh;
	int nBitsOffset = sizeof(BITMAPFILEHEADER) + BMIH.biSize;
	LONG lImageSize = BMIH.biSizeImage;
	LONG lFileSize = nBitsOffset + lImageSize;
	bmfh.bfType = 'B' + ('M' << 8);
	bmfh.bfOffBits = nBitsOffset;
	bmfh.bfSize = lFileSize;
	bmfh.bfReserved1 = bmfh.bfReserved2 = 0;

	//Write the bitmap file header
	UINT nWrittenFileHeaderSize = fwrite(&bmfh, 1,
		sizeof(BITMAPFILEHEADER), pFile);

	//And then the bitmap info header
	UINT nWrittenInfoHeaderSize = fwrite(&BMIH,
		1, sizeof(BITMAPINFOHEADER), pFile);

	//Finally, write the image data itself 

	//-- the data represents our drawing

	UINT nWrittenDIBDataSize =
		fwrite(pixels, 1, lImageSize, pFile);

	fclose(pFile);

	delete pixels;

	return true;
}

int SaveTGA(string filename) {
	// we will store the image data here
	unsigned char *pixels;
	// the thingy we use to write files
	FILE * shot;
	// we get the width/height of the screen into this array
	int screenStats[4];

	// get the width/height of the window
	glGetIntegerv(GL_VIEWPORT, screenStats);

	// generate an array large enough to hold the pixel data 
	// (width*height*bytesPerPixel)
	pixels = new unsigned char[screenStats[2] * screenStats[3] * 3];
	// read in the pixel data, TGA's pixels are BGR aligned
	glReadPixels(0, 0, screenStats[2], screenStats[3], GL_BGR,
		GL_UNSIGNED_BYTE, pixels);

	// open the file for writing. If unsucessful, return 1
	//if ((fopen_s(&shot, "shot.tga", "wb")) == NULL) return 1;
	//std::ostringstream STEPNUM; STEPNUM << gStepNum;

	//string filename = "imgs\\imgdata" + to_string(gStepNum) + ".tga";
	//string address = OUTPUT_DATA + filename;
	fopen_s(&shot, filename.c_str(), "wb");

	// this is the tga header it must be in the beginning of 
	// every (uncompressed) .tga
	unsigned char TGAheader[12] = { 0,0,2,0,0,0,0,0,0,0,0,0 };
	// the header that is used to get the dimensions of the .tga
	// header[1]*256+header[0] - width
	// header[3]*256+header[2] - height
	// header[4] - bits per pixel
	// header[5] - ?
	unsigned char header[6] = { ((int)(screenStats[2] % 256)),
		((int)(screenStats[2] / 256)),
		((int)(screenStats[3] % 256)),
		((int)(screenStats[3] / 256)),24,0 };

	// write out the TGA header
	fwrite(TGAheader, sizeof(unsigned char), 12, shot);
	// write out the header
	fwrite(header, sizeof(unsigned char), 6, shot);
	// write the pixels
	fwrite(pixels, sizeof(unsigned char),
		screenStats[2] * screenStats[3] * 3, shot);

	// close the file
	fclose(shot);
	// free the memory
	delete[] pixels;

	// return success
	return 0;
}