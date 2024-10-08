//
// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) Contributors to the OpenEXR Project.
//

#ifdef NDEBUG
#    undef NDEBUG
#endif

#include "compareB44.h"
#include "compareDwa.h"

#include <IlmThread.h>
#include <ImfArray.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <ImfHeader.h>
#include <ImfMultiPartOutputFile.h>
#include <ImfOutputFile.h>
#include <ImfOutputPart.h>
#include <ImfPartType.h>
#include <ImfCRgbaFile.h>
#include <ImfThreading.h>

#include <string>
#include <vector>

#include <assert.h>
#include <stdio.h>

using namespace OPENEXR_IMF_NAMESPACE; // should this be OPENEXR_IMF_INTERNAL_NAMESPACE?
using namespace std;
using namespace IMATH_NAMESPACE;

namespace
{

void
imfRgbaMethods ()
{
    //
    // Verify that the constructors and the assignment
    // operator for struct ImfRgba work.
    //

    ImfRgba defaultImfRgba; // Test default initialization
    assert (defaultImfRgba.r == 0 && defaultImfRgba.g == 0 && 
                defaultImfRgba.b == 0 && defaultImfRgba.a == 0);

    // TO DO : The values must be unsigned shorts, are they?
    ImfRgba x = {2, 3, 4};
    assert (x.r == 2 && x.g == 3 && x.b == 4 && x.a == 0);

    ImfRgba y = {5, 6, 7, 0};
    assert (y.r == 5 && y.g == 6 && y.b == 7 && y.a == 0);

    ImfRgba z; 

    z = x;
    assert (z.r == 2 && z.g == 3 && z.b == 4 && z.a == 0);

    z = y;
    assert (z.r == 5 && z.g == 6 && z.b == 7 && z.a == 0);

    // ImfRgba w (z);
    // assert (w.r == 5.f && w.g == 6.f && w.b == 7.f && w.a == 0.f);
}

void
fillPixels (ImfRgba* pixels, int w, int h)
{
    // for (int y = 0; y < h; ++y)
    // {
    //     for (int x = 0; x < w; ++x)
    //     {
    //         // ImfRgba& p = pixels[y][x];
    //         pixels[y][x].r = 0.5 + 0.5 * sin (0.1 * x + 0.1 * y);
    //         pixels[y][x].g = 0.5 + 0.5 * sin (0.1 * x + 0.2 * y);
    //         pixels[y][x].b = 0.5 + 0.5 * sin (0.1 * x + 0.3 * y);
    //         pixels[y][x].a = (pixels[y][x].r + pixels[y][x].b + pixels[y][x].g) / 3.0;
    //     }
    // }
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int index = y * w + x; // Calculate 1D index
            pixels[index].r = 0.5 + 0.5 * sin(0.1 * x + 0.1 * y);
            pixels[index].g = 0.5 + 0.5 * sin(0.1 * x + 0.2 * y);
            pixels[index].b = 0.5 + 0.5 * sin(0.1 * x + 0.3 * y);
            pixels[index].a = (pixels[index].r + pixels[index].b + pixels[index].g) / 3.0;
        }
    }
}

void convertImfRgbaToRgba(const Array2D<ImfRgba>& ImfRgbaPixels, Array2D<Rgba>& RgbaPixels, int w, int h)
{
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            RgbaPixels[y][x].r = ImfRgbaPixels[y][x].r;
            RgbaPixels[y][x].g = ImfRgbaPixels[y][x].g;
            RgbaPixels[y][x].b = ImfRgbaPixels[y][x].b;
            RgbaPixels[y][x].a = ImfRgbaPixels[y][x].a;
        }
    }
}

void
writeReadCRGBA (
    const char           fileName[],
    int                  width,
    int                  height,
    const ImfRgba* p1,
    RgbaChannels         channels,
    LineOrder            lorder,
    Compression          comp)
{
    //
    // Save the selected channels of RGBA image p1; save the
    // scan lines in the specified order.  Read the image back
    // from the file, and compare the data with the original.
    //

    cout << "channels " << ((channels & IMF_WRITE_R) ? "R" : "")
         << ((channels & IMF_WRITE_G) ? "G" : "")
         << ((channels & IMF_WRITE_B) ? "B" : "")
         << ((channels & IMF_WRITE_A) ? "A" : "") << ", line order " << lorder
         << ", compression " << comp << endl;

    // Header header (width, height);
    // header.lineOrder ()   = lorder;
    // header.compression () = comp;

    ImfHeader* headerPtr = ImfNewHeader(); // TO DO : delete headerPtr?
    ImfHeaderSetLineOrder(headerPtr, lorder);
    ImfHeaderSetCompression(headerPtr, comp);
    ImfHeaderSetDataWindow(headerPtr, 0, 0, width - 1, height - 1); // TO DO : Should it be width - 1, height - 1?
    ImfHeaderSetDisplayWindow(headerPtr, 0, 0, width - 1, height - 1);
    // TO DO : set dataWindow using a setBox2iAttribute function? or is Box2i meant for generic attributes?
    // TO DO : Delete headerPtr?

    // Get dataWindow
    int xMinDataWindowHeader, yMinDataWindowHeader, xMaxDataWindowHeader, yMaxDataWindowHeader;
    ImfHeaderDataWindow(headerPtr, &xMinDataWindowHeader, &yMinDataWindowHeader, &xMaxDataWindowHeader, &yMaxDataWindowHeader);

    // Get displayWindow
    int xMinDisplayWindowHeader, yMinDisplayWindowHeader, xMaxDisplayWindowHeader, yMaxDisplayWindowHeader;
    ImfHeaderDisplayWindow(headerPtr, &xMinDisplayWindowHeader, &yMinDisplayWindowHeader, &xMaxDisplayWindowHeader, &yMaxDisplayWindowHeader);

    // Get screenWindowCenter
    float xScreenWindowCenterHeader, yScreenWindowCenterHeader;
    ImfHeaderScreenWindowCenter(headerPtr, &xScreenWindowCenterHeader, &yScreenWindowCenterHeader);

    cout << "writing ";
    cout.flush ();

    {
        // remove (fileName);
        // RgbaOutputFile out (fileName, header, channels);
        // out.setFrameBuffer (&p1[0][0], 1, width); 
        // out.writePixels (height);

        ImfOutputFile* out = ImfOpenOutputFile(fileName, headerPtr, channels);
        ImfOutputSetFrameBuffer(out, (ImfRgba*) p1, 1, width); // TO DO : Is this correct?
        ImfOutputWritePixels(out, height); // TO DO : Is this correct?
        int outputFileClosed = ImfCloseOutputFile(out); // TO DO : delete file ptr correct?
        assert (outputFileClosed == 1);
    }


    cout << "reading ";
    cout.flush ();

    // {
    //     // RgbaInputFile in (fileName);
    //     ImfInputFile* in = ImfOpenInputFile(fileName); // TO DO : delete ptr?
    //     // const Box2i&  dw = in.dataWindow ();
    //     const ImfHeader* inputFileHeaderPtr = ImfInputHeader(in); // TO DO : Where does this get deleted?
    //     int xMin, yMin, xMax, yMax;
    //     ImfHeaderDataWindow(inputFileHeaderPtr, &xMin, &yMin, &xMax, &yMax);

    //     int w  = xMax - xMin + 1; // dw.max.x - dw.min.x + 1;
    //     int h  = yMax - yMin + 1;// dw.max.y - dw.min.y + 1;
    //     int dx = xMin; // dw.min.x;
    //     int dy = yMin; // dw.min.y;

    //     // Array2D<Rgba> p2 (h, w);
    //     Array2D<ImfRgba> p2(h, w);
    //     for (int y = 0; y < h; ++y)
    //     {
    //         for (int x = 0; x < w; ++x)
    //         {
    //             p2[y][x] = {0, 0, 0, 0};
    //             // p2[y][x].r = 0;
    //             // p2[y][x].g = 0;
    //             // p2[y][x].b = 0;
    //             // p2[y][x].a = 0;
    //         }
    //     }
    //     // in.setFrameBuffer (&p2[-dy][-dx], 1, w);
    //     int successSetFrameBuffer = ImfInputSetFrameBuffer(in, &p2[-dy][-dx], 1, w); // TO DO : Check if success == 1?
    //     assert(successSetFrameBuffer == 1);
    //     // in.readPixels (dw.min.y, dw.max.y);
    //     int successReadPixels = ImfInputReadPixels(in, yMin, yMax); // TO DO : Check if success == 1?
    //     assert(successReadPixels == 1);

    //     // TO DO : Pending checks for display window, data window, screen window center (DONE)

    //     // Get inputFile dataWindow
    //     int xMinDataWindowInputFileHeader, yMinDataWindowInputFileHeader, xMaxDataWindowInputFileHeader, yMaxDataWindowInputFileHeader;
    //     ImfHeaderDataWindow(headerPtr, &xMinDataWindowInputFileHeader, &yMinDataWindowInputFileHeader, &xMaxDataWindowInputFileHeader, &yMaxDataWindowInputFileHeader);

    //     // Get inputFile displayWindow
    //     int xMinDisplayWindowInputFileHeader, yMinDisplayWindowInputFileHeader, xMaxDisplayWindowInputFileHeader, yMaxDisplayWindowInputFileHeader;
    //     ImfHeaderDisplayWindow(headerPtr, &xMinDisplayWindowInputFileHeader, &yMinDisplayWindowInputFileHeader, &xMaxDisplayWindowInputFileHeader, &yMaxDisplayWindowInputFileHeader);

    //     // Get inputFile screenWindowCenter
    //     float xScreenWindowCenterInputFileHeader, yScreenWindowCenterInputFileHeader;
    //     ImfHeaderScreenWindowCenter(headerPtr, &xScreenWindowCenterInputFileHeader, &yScreenWindowCenterInputFileHeader);

        
    //     // assert (in.displayWindow () == header.displayWindow ());
    //     assert ((xMinDisplayWindowHeader == xMinDisplayWindowInputFileHeader) && (yMinDisplayWindowHeader == yMinDisplayWindowInputFileHeader) && (xMaxDisplayWindowHeader == xMaxDisplayWindowInputFileHeader) && (yMaxDisplayWindowHeader == yMaxDisplayWindowInputFileHeader));
    //     // assert (in.dataWindow () == header.dataWindow ());
    //     assert ((xMinDataWindowHeader == xMinDataWindowInputFileHeader) && (yMinDataWindowHeader == yMinDataWindowInputFileHeader) && (xMaxDataWindowHeader == xMaxDataWindowInputFileHeader) && (yMaxDataWindowHeader == yMaxDataWindowInputFileHeader));
    //     // assert (in.pixelAspectRatio () == header.pixelAspectRatio ());
    //     assert (ImfHeaderPixelAspectRatio(inputFileHeaderPtr) == ImfHeaderPixelAspectRatio(headerPtr));
    //     // assert (in.screenWindowCenter () == header.screenWindowCenter ());
    //     assert ((xScreenWindowCenterHeader == xScreenWindowCenterInputFileHeader) && (yScreenWindowCenterHeader == yScreenWindowCenterInputFileHeader));
    //     // assert (in.screenWindowWidth () == header.screenWindowWidth ());
    //     assert (ImfHeaderScreenWindowWidth(inputFileHeaderPtr) == ImfHeaderScreenWindowWidth(headerPtr));
    //     // assert (in.lineOrder () == header.lineOrder ());
    //     assert (ImfHeaderLineOrder(inputFileHeaderPtr) == ImfHeaderLineOrder(headerPtr));
    //     // assert (in.compression () == header.compression ())
    //     assert (ImfHeaderCompression(inputFileHeaderPtr) == ImfHeaderCompression(headerPtr));
    //     // assert (in.channels () == channels);
    //     assert(ImfInputChannels(in) == channels);

    //     // if (in.compression () == B44_COMPRESSION ||
    //     //     in.compression () == B44A_COMPRESSION)
    //     // {
    //     //     compareB44 (width, height, p1, p2, channels);
    //     // }
    //     // else if (
    //     //     in.compression () == DWAA_COMPRESSION ||
    //     //     in.compression () == DWAB_COMPRESSION)
    //     // {
    //     //     compareDwa (width, height, p1, p2, channels);
    //     // }
    //     // else
    //     // {
    //     //     for (int y = 0; y < h; ++y)
    //     //     {
    //     //         for (int x = 0; x < w; ++x)
    //     //         {
    //     //             if (channels & WRITE_R)
    //     //                 assert (p2[y][x].r == p1[y][x].r);
    //     //             else
    //     //                 assert (p2[y][x].r == 0);

    //     //             if (channels & WRITE_G)
    //     //                 assert (p2[y][x].g == p1[y][x].g);
    //     //             else
    //     //                 assert (p2[y][x].g == 0);

    //     //             if (channels & WRITE_B)
    //     //                 assert (p2[y][x].b == p1[y][x].b);
    //     //             else
    //     //                 assert (p2[y][x].b == 0);

    //     //             if (channels & WRITE_A)
    //     //                 assert (p2[y][x].a == p1[y][x].a);
    //     //             else
    //     //                 assert (p2[y][x].a == 1);
    //     //         }
    //     //     }
    //     // }

    //     if (ImfHeaderCompression(inputFileHeaderPtr) == IMF_B44_COMPRESSION ||
    //         ImfHeaderCompression(inputFileHeaderPtr) == IMF_B44A_COMPRESSION)
    //     {
    //         // Array2D<Rgba> RgbaPixelsp1(h, w);
    //         // convertImfRgbaToRgba(p1, RgbaPixelsp1, w, h);
    //         // Array2D<Rgba> RgbaPixelsp2(h, w);
    //         // convertImfRgbaToRgba(p2, RgbaPixelsp2, w, h);
    //         // compareB44(width, height, RgbaPixelsp1, RgbaPixelsp2, channels); // does this need to take a reference?
    //     }
    //     else if (
    //         ImfHeaderCompression(inputFileHeaderPtr) == IMF_DWAA_COMPRESSION ||
    //         ImfHeaderCompression(inputFileHeaderPtr) == IMF_DWAB_COMPRESSION)
    //     {
    //         // Array2D<Rgba> RgbaPixelsp1(h, w);
    //         // convertImfRgbaToRgba(p1, RgbaPixelsp1, w, h);
    //         // Array2D<Rgba> RgbaPixelsp2(h, w);
    //         // convertImfRgbaToRgba(p2, RgbaPixelsp2, w, h);
    //         // compareDwa(width, height, RgbaPixelsp1, RgbaPixelsp2, channels);
    //     }
    //     else
    //     {
    //         for (int y = 0; y < h; ++y)
    //         {
    //             for (int x = 0; x < w; ++x)
    //             {
    //                 cout << "p1[y][x].rgba = (" << p1[y][x].r << ", " << 
    //                 p1[y][x].g << ", " << p1[y][x].b << ", " << p1[y][x].a << endl;
    //                 cout << "p2[y][x].rgba = (" << p2[y][x].r << ", " << 
    //                 p2[y][x].g << ", " << p2[y][x].b << ", " << p2[y][x].a << endl;
                    
    //                 if (channels & IMF_WRITE_R)
    //                     assert (p2[y][x].r == p1[y][x].r);
    //                 else
    //                     assert (p2[y][x].r == 0);

    //                 if (channels & IMF_WRITE_G)
    //                     assert (p2[y][x].g == p1[y][x].g);
    //                 else
    //                     assert (p2[y][x].g == 0);

    //                 if (channels & IMF_WRITE_B)
    //                     assert (p2[y][x].b == p1[y][x].b);
    //                 else
    //                     assert (p2[y][x].b == 0);

    //                 if (channels & IMF_WRITE_A)
    //                     assert (p2[y][x].a == p1[y][x].a);
    //                 else
    //                     assert (p2[y][x].a == 0);
    //             }
    //         }
    //     }

    //     // ImfDeleteHeader(inputFileHeaderPtr);
    //     int inputFileClosed = ImfCloseInputFile(in);
    //     assert (inputFileClosed == 1);
    // }
    ImfDeleteHeader(headerPtr); // TO DO : should these be deleted earlier?
    // remove (fileName); // TO DO : This is not required? But I dunno what this function does
}

} // namespace

void
testCRgba (const std::string& tempDir)
{
    std::string tempDir1 = "/Users/megs/Documents/aswf/take2_copy/";
    try
    {
        std::cout << "Testing the RGBA image interface....." << tempDir1 + "imf_test_rgba.exr" << std::endl;

        imfRgbaMethods ();

        const int W = 237;
        const int H = 119;

        // Array2D<ImfRgba> p1 (H, W);
        ImfRgba* p1 = (ImfRgba*) malloc(H * W * sizeof(ImfRgba));
        // fillPixels ((ImfRgba*)p1, W, H);
        for (int i=0;i<H;i++) {
            for (int j=0;j<W;j++) {
                ImfFloatToHalf (0.5, &p1[i * W + j].r); // 3;//0.5 + 0.5 * sin(0.1 * j + 0.1 * i);
                ImfFloatToHalf (0.5, &p1[i * W + j].g);
                ImfFloatToHalf (0.5, &p1[i * W + j].b);
                ImfFloatToHalf (0.5, &p1[i * W + j].a);
                // p1[i * W + j].r = 1;//0.5 + 0.5 * sin(0.1 * j + 0.1 * i);
                // p1[i * W + j].g = 1;//0.5 + 0.5 * sin(0.1 * j + 0.2 * i);
                // p1[i * W + j].b = 1;//0.5 + 0.5 * sin(0.1 * j + 0.3 * i);
                // p1[i * W + j].a = 0;
                cout << "p1[" << i << "][" << j << "].rgba = (" << p1[i * W + j].r << ", " << 
                p1[i * W + j].g << ", " << p1[i * W + j].b << ", " << p1[i * W + j].a << ")" << endl;
            }
        }

        int maxThreads = ILMTHREAD_NAMESPACE::supportsThreads () ? 3 : 0;

        for (int n = 0; n <= maxThreads; ++n) // Vary num threads 0 to max(3)
        {
            if (ILMTHREAD_NAMESPACE::supportsThreads ())
            {
                setGlobalThreadCount (n);
                cout << "\nnumber of threads: " << globalThreadCount () << endl;
            }

            for (int lorder = 0; lorder < RANDOM_Y; ++lorder) // Vary line order
            {
                for (int comp = 0; comp < NUM_COMPRESSION_METHODS; ++comp) // Vary compression
                {
                    writeReadCRGBA (
                        (tempDir1 + "imf_test_crgba.exr").c_str (),
                        W,
                        H,
                        (ImfRgba*) p1,
                        WRITE_RGBA,
                        LineOrder (lorder),
                        Compression (comp));

                    // writeReadCRGBA (
                    //     (tempDir + "imf_test_rgba.exr").c_str (),
                    //     W,
                    //     H,
                    //     p1,
                    //     WRITE_RGB,
                    //     LineOrder (lorder),
                    //     Compression (comp));

                    // writeReadCRGBA (
                    //     (tempDir + "imf_test_rgba.exr").c_str (),
                    //     W,
                    //     H,
                    //     p1,
                    //     WRITE_A,
                    //     LineOrder (lorder),
                    //     Compression (comp));

                    // writeReadCRGBA (
                    //     (tempDir + "imf_test_rgba.exr").c_str (),
                    //     W,
                    //     H,
                    //     p1,
                    //     RgbaChannels (WRITE_R | WRITE_B),
                    //     LineOrder (lorder),
                    //     Compression (comp));
                }
            }

            // writeReadIncomplete (tempDir);
        }

        // writeReadLayers (tempDir, false);
        // writeReadLayers (tempDir, true);

        cout << "ok\n" << endl;
        free(p1);
    }
    catch (const std::exception& e)
    {
        cerr << "ERROR -- caught exception: " << e.what () << endl;
        assert (false);
    }
}
