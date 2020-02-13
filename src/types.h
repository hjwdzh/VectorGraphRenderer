#ifndef VECTORGRAPH_RENDERER_TYPES_H_
#define VECTORGRAPH_RENDERER_TYPES_H_

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

struct FaceIDPair
{
	long long operator()( const long long& lhs, const long long& rhs ) const {
		return lhs * ((long long)(1 << 16)) * ((long long)(1 << 16)) + rhs;
	}
};

typedef CGAL::Exact_rational									K;
typedef CGAL::Cartesian<K>                   					Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::X_monotone_curve_2                            Segment_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, long long>            Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>                     Arrangement_2;
typedef CGAL::Arr_face_overlay_traits<Arrangement_2,
                                      Arrangement_2,
                                      Arrangement_2,
                                      FaceIDPair >  Overlay_traits;


#endif