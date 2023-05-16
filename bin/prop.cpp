#include "global.hpp"
#include "aliases.hpp"
#include "operations.hpp"

#define EXTERN_PROP
 #include "prop.hpp"

namespace qprop
{
    void set_ins()
    {
        if(ntypes==1) ins_list={LO};
        nins=ins_list.size();
    }
}
namespace jprop
{
    void set_ins()
    {
        if(ntypes==1) ins_list={LO};
        nins=ins_list.size();
    }
}
namespace lprop
{
    void set_ins()
    {
        if(ntypes==1) ins_list={LO};
        else ins_list={LO,F};
        nins=ins_list.size();
    }
}

void oper_t::build_prop(const vvvprop_t &prop, vvvprop_t &jprop)
{
#pragma omp parallel for collapse(2)
    for(int mr=0;mr<nmr;mr++)
        for(int ijack=0;ijack<njacks;ijack++)
        {
            if(ntypes==1)
            {
                jprop[jprop::LO ][ijack][mr] += prop[ijack][qprop::LO ][mr]; // Leading order
            }   
        }
}


// invert the propagator
vvprop_t invert_jprop(const vvprop_t &jprop)
{
    vvprop_t jprop_inv(valarray<prop_t>(prop_t::Zero(),nmr),njacks);
    
#pragma omp parallel for collapse(2)
    for(int ijack=0;ijack<njacks;ijack++)
        for(int mr=0;mr<nmr;mr++)
            jprop_inv[ijack][mr]=jprop[ijack][mr].inverse();
    
    return jprop_inv;
}