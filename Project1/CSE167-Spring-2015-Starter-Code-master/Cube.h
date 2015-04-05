#ifndef CSE167_Cube_h
#define CSE167_Cube_h

#include "Drawable.h"

class Cube : public Drawable
{
    
public:
    
    float size;
    
    Cube(float);
    virtual ~Cube(void);
    
    virtual void draw(DrawData&);
    virtual void update(UpdateData&);
    
    void spin(float); 
    void translate(float, float, float);
    void scale(float);
    void reset();
    void orbit(float);
    inline void setToggle(bool value){toggle = value;}
    inline bool getToggle(){return toggle;}
private:
	bool toggle = false;
};

#endif

