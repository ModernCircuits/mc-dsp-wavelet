#pragma once

#include <mc/core/array.hpp>

namespace mc::dsp {

template<typename T>
inline constexpr auto const h1 = Array<T, 10>{
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.70710678118654752440084436210),
    static_cast<T>(0.70710678118654752440084436210),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm111 = Array<T, 2>{
    static_cast<T>(0.70710678118654752440084436210),
    static_cast<T>(0.70710678118654752440084436210),
};

template<typename T>
inline constexpr auto const hm113 = Array<T, 6>{
    static_cast<T>(-0.0883883476483184405501055452631),
    static_cast<T>(0.0883883476483184405501055452631),
    static_cast<T>(0.70710678118654752440084436210),
    static_cast<T>(0.70710678118654752440084436210),
    static_cast<T>(0.0883883476483184405501055452631),
    static_cast<T>(-0.0883883476483184405501055452631),
};

template<typename T>
inline constexpr auto const hm115 = Array<T, 10>{
    static_cast<T>(0.0165728151840597076031447897368),
    static_cast<T>(-0.0165728151840597076031447897368),
    static_cast<T>(-0.1215339780164378557563951247368),
    static_cast<T>(0.1215339780164378557563951247368),
    static_cast<T>(0.70710678118654752440084436210),
    static_cast<T>(0.70710678118654752440084436210),
    static_cast<T>(0.1215339780164378557563951247368),
    static_cast<T>(-0.1215339780164378557563951247368),
    static_cast<T>(-0.0165728151840597076031447897368),
    static_cast<T>(0.0165728151840597076031447897368),
};

template<typename T>
inline constexpr auto const h2 = Array<T, 18>{
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.3535533905932737622004221810524),
    static_cast<T>(0.7071067811865475244008443621048),
    static_cast<T>(0.3535533905932737622004221810524),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm222 = Array<T, 6>{
    static_cast<T>(-0.1767766952966368811002110905262),
    static_cast<T>(0.3535533905932737622004221810524),
    static_cast<T>(1.0606601717798212866012665431573),
    static_cast<T>(0.3535533905932737622004221810524),
    static_cast<T>(-0.1767766952966368811002110905262),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm224 = Array<T, 10>{
    static_cast<T>(0.0331456303681194152062895794737),
    static_cast<T>(-0.0662912607362388304125791589473),
    static_cast<T>(-0.1767766952966368811002110905262),
    static_cast<T>(0.4198446513295125926130013399998),
    static_cast<T>(0.9943689110435824561886873842099),
    static_cast<T>(0.4198446513295125926130013399998),
    static_cast<T>(-0.1767766952966368811002110905262),
    static_cast<T>(-0.0662912607362388304125791589473),
    static_cast<T>(0.0331456303681194152062895794737),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm226 = Array<T, 14>{
    static_cast<T>(-0.0069053396600248781679769957237),
    static_cast<T>(0.0138106793200497563359539914474),
    static_cast<T>(0.0469563096881691715422435709210),
    static_cast<T>(-0.1077232986963880994204411332894),
    static_cast<T>(-0.1698713556366120029322340948025),
    static_cast<T>(0.4474660099696121052849093228945),
    static_cast<T>(0.9667475524034829435167794013152),
    static_cast<T>(0.4474660099696121052849093228945),
    static_cast<T>(-0.1698713556366120029322340948025),
    static_cast<T>(-0.1077232986963880994204411332894),
    static_cast<T>(0.0469563096881691715422435709210),
    static_cast<T>(0.0138106793200497563359539914474),
    static_cast<T>(-0.0069053396600248781679769957237),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm228 = Array<T, 18>{
    static_cast<T>(0.0015105430506304420992449678146),
    static_cast<T>(-0.0030210861012608841984899356291),
    static_cast<T>(-0.0129475118625466465649568669819),
    static_cast<T>(0.0289161098263541773284036695929),
    static_cast<T>(0.0529984818906909399392234421792),
    static_cast<T>(-0.1349130736077360572068505539514),
    static_cast<T>(-0.1638291834340902345352542235443),
    static_cast<T>(0.4625714404759165262773590010400),
    static_cast<T>(0.9516421218971785225243297231697),
    static_cast<T>(0.4625714404759165262773590010400),
    static_cast<T>(-0.1638291834340902345352542235443),
    static_cast<T>(-0.1349130736077360572068505539514),
    static_cast<T>(0.0529984818906909399392234421792),
    static_cast<T>(0.0289161098263541773284036695929),
    static_cast<T>(-0.0129475118625466465649568669819),
    static_cast<T>(-0.0030210861012608841984899356291),
    static_cast<T>(0.0015105430506304420992449678146),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const h3 = Array<T, 20>{
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.1767766952966368811002110905262),
    static_cast<T>(0.5303300858899106433006332715786),
    static_cast<T>(0.5303300858899106433006332715786),
    static_cast<T>(0.1767766952966368811002110905262),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm331 = Array<T, 4>{
    static_cast<T>(-0.3535533905932737622004221810524),
    static_cast<T>(1.0606601717798212866012665431573),
    static_cast<T>(1.0606601717798212866012665431573),
    static_cast<T>(-0.3535533905932737622004221810524),
};

template<typename T>
inline constexpr auto const hm333 = Array<T, 8>{
    static_cast<T>(0.0662912607362388304125791589473),
    static_cast<T>(-0.1988737822087164912377374768420),
    static_cast<T>(-0.1546796083845572709626847042104),
    static_cast<T>(0.9943689110435824561886873842099),
    static_cast<T>(0.9943689110435824561886873842099),
    static_cast<T>(-0.1546796083845572709626847042104),
    static_cast<T>(-0.1988737822087164912377374768420),
    static_cast<T>(0.0662912607362388304125791589473),
};

template<typename T>
inline constexpr auto const hm335 = Array<T, 12>{
    static_cast<T>(-0.0138106793200497563359539914474),
    static_cast<T>(0.0414320379601492690078619743421),
    static_cast<T>(0.0524805814161890740766251675000),
    static_cast<T>(-0.2679271788089652729175074340788),
    static_cast<T>(-0.0718155324642587329469607555263),
    static_cast<T>(0.9667475524034829435167794013152),
    static_cast<T>(0.9667475524034829435167794013152),
    static_cast<T>(-0.0718155324642587329469607555263),
    static_cast<T>(-0.2679271788089652729175074340788),
    static_cast<T>(0.0524805814161890740766251675000),
    static_cast<T>(0.0414320379601492690078619743421),
    static_cast<T>(-0.0138106793200497563359539914474),
};

template<typename T>
inline constexpr auto const hm337 = Array<T, 16>{
    static_cast<T>(0.0030210861012608841984899356291),
    static_cast<T>(-0.0090632583037826525954698068873),
    static_cast<T>(-0.0168317654213106405344439270765),
    static_cast<T>(0.0746639850740189951912512662623),
    static_cast<T>(0.0313329787073628846871956180962),
    static_cast<T>(-0.3011591259228349991008967259990),
    static_cast<T>(-0.0264992409453454699696117210896),
    static_cast<T>(0.9516421218971785225243297231697),
    static_cast<T>(0.9516421218971785225243297231697),
    static_cast<T>(-0.0264992409453454699696117210896),
    static_cast<T>(-0.3011591259228349991008967259990),
    static_cast<T>(0.0313329787073628846871956180962),
    static_cast<T>(0.0746639850740189951912512662623),
    static_cast<T>(-0.0168317654213106405344439270765),
    static_cast<T>(-0.0090632583037826525954698068873),
    static_cast<T>(0.0030210861012608841984899356291),
};

template<typename T>
inline constexpr auto const hm339 = Array<T, 20>{
    static_cast<T>(-0.0006797443727836989446602355165),
    static_cast<T>(0.0020392331183510968339807065496),
    static_cast<T>(0.0050603192196119810324706421788),
    static_cast<T>(-0.0206189126411055346546938106687),
    static_cast<T>(-0.0141127879301758447558029850103),
    static_cast<T>(0.0991347824942321571990197448581),
    static_cast<T>(0.0123001362694193142367090236328),
    static_cast<T>(-0.3201919683607785695513833204624),
    static_cast<T>(0.0020500227115698857061181706055),
    static_cast<T>(0.9421257006782067372990864259380),
    static_cast<T>(0.9421257006782067372990864259380),
    static_cast<T>(0.0020500227115698857061181706055),
    static_cast<T>(-0.3201919683607785695513833204624),
    static_cast<T>(0.0123001362694193142367090236328),
    static_cast<T>(0.0991347824942321571990197448581),
    static_cast<T>(-0.0141127879301758447558029850103),
    static_cast<T>(-0.0206189126411055346546938106687),
    static_cast<T>(0.0050603192196119810324706421788),
    static_cast<T>(0.0020392331183510968339807065496),
    static_cast<T>(-0.0006797443727836989446602355165),
};

template<typename T>
inline constexpr auto const h4 = Array<T, 10>{
    static_cast<T>(0.0),
    static_cast<T>(-0.064538882628697058),
    static_cast<T>(-0.040689417609164058),
    static_cast<T>(0.41809227322161724),
    static_cast<T>(0.7884856164055829),
    static_cast<T>(0.41809227322161724),
    static_cast<T>(-0.040689417609164058),
    static_cast<T>(-0.064538882628697058),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm444 = Array<T, 10>{
    static_cast<T>(0.03782845550726404),
    static_cast<T>(-0.023849465019556843),
    static_cast<T>(-0.11062440441843718),
    static_cast<T>(0.37740285561283066),
    static_cast<T>(0.85269867900889385),
    static_cast<T>(0.37740285561283066),
    static_cast<T>(-0.11062440441843718),
    static_cast<T>(-0.023849465019556843),
    static_cast<T>(0.03782845550726404),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const h5 = Array<T, 12>{
    static_cast<T>(0.013456709459118716),
    static_cast<T>(-0.0026949668801115071),
    static_cast<T>(-0.13670658466432914),
    static_cast<T>(-0.093504697400938863),
    static_cast<T>(0.47680326579848425),
    static_cast<T>(0.89950610974864842),
    static_cast<T>(0.47680326579848425),
    static_cast<T>(-0.093504697400938863),
    static_cast<T>(-0.13670658466432914),
    static_cast<T>(-0.0026949668801115071),
    static_cast<T>(0.013456709459118716),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm555 = Array<T, 12>{
    static_cast<T>(0.0),
    static_cast<T>(0.03968708834740544),
    static_cast<T>(0.0079481086372403219),
    static_cast<T>(-0.054463788468236907),
    static_cast<T>(0.34560528195603346),
    static_cast<T>(0.73666018142821055),
    static_cast<T>(0.34560528195603346),
    static_cast<T>(-0.054463788468236907),
    static_cast<T>(0.0079481086372403219),
    static_cast<T>(0.03968708834740544),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const h6 = Array<T, 18>{
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.014426282505624435),
    static_cast<T>(0.014467504896790148),
    static_cast<T>(-0.078722001062628819),
    static_cast<T>(-0.040367979030339923),
    static_cast<T>(0.41784910915027457),
    static_cast<T>(0.75890772945365415),
    static_cast<T>(0.41784910915027457),
    static_cast<T>(-0.040367979030339923),
    static_cast<T>(-0.078722001062628819),
    static_cast<T>(0.014467504896790148),
    static_cast<T>(0.014426282505624435),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
    static_cast<T>(0.0),
};

template<typename T>
inline constexpr auto const hm668 = Array<T, 18>{
    static_cast<T>(0.0019088317364812906),
    static_cast<T>(-0.0019142861290887667),
    static_cast<T>(-0.016990639867602342),
    static_cast<T>(0.01193456527972926),
    static_cast<T>(0.04973290349094079),
    static_cast<T>(-0.077263173167204144),
    static_cast<T>(-0.09405920349573646),
    static_cast<T>(0.42079628460982682),
    static_cast<T>(0.82592299745840225),
    static_cast<T>(0.42079628460982682),
    static_cast<T>(-0.09405920349573646),
    static_cast<T>(-0.077263173167204144),
    static_cast<T>(0.04973290349094079),
    static_cast<T>(0.01193456527972926),
    static_cast<T>(-0.016990639867602342),
    static_cast<T>(-0.0019142861290887667),
    static_cast<T>(0.0019088317364812906),
    static_cast<T>(0.0),
};
}  // namespace mc::dsp
