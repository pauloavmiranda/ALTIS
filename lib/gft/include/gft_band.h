#ifndef _GFT_BAND_H_
#define _GFT_BAND_H_

#include "gft_common.h"
#include "gft_image32.h"

namespace gft{

  namespace Band{

    /**
     * \brief Creates a circular template for band shape constraints.
     *
     * @param size The dimensions (length and width) of the template in pixels, which must be an odd number. 
     */
    gft::sImage32 *CircularTemplate(int size);

    /**
     * \brief Creates a rectangular template for band shape constraints.
     *
     * @param size The dimensions (length and width) of the template in pixels, which must be an odd number. 
     * @param ratio The aspect ratio width/height.
     */
    gft::sImage32 *RectangularTemplate(int size, float ratio);

    /**
     * \brief Creates an elliptical template for band shape constraints.
     *
     * @param size The dimensions (length and width) of the template in pixels, which must be an odd number. 
     * @param ratio The aspect ratio width/height.
     */
    gft::sImage32 *EllipticalTemplate(int size, float ratio);

    /**
     * \brief Creates a template by Gielis equation for band shape constraints.
     *
     * @param size The dimensions (length and width) of the template in pixels, which must be an odd number. 
     */
    gft::sImage32 *GielisEquationTemplate(int size,
					  float A,
					  float B,
					  float n1,
					  float n2,
					  float n3,
					  int m);
    
  } //end Band namespace

} //end gft namespace



#endif
