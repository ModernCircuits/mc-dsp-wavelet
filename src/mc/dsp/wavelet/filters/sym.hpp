// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/array.hpp>

namespace mc {

template<typename T>
inline constexpr auto sym2 = Array<T, 4>{
    static_cast<T>(0.48296291314469025),
    static_cast<T>(0.83651630373746899),
    static_cast<T>(0.22414386804185735),
    static_cast<T>(-0.12940952255092145),
};

template<typename T>
inline constexpr auto sym3 = Array<T, 6>{
    static_cast<T>(0.33267055295095688),
    static_cast<T>(0.80689150931333875),
    static_cast<T>(0.45987750211933132),
    static_cast<T>(-0.13501102001039084),
    static_cast<T>(-0.085441273882241486),
    static_cast<T>(0.035226291882100656),
};

template<typename T>
inline constexpr auto sym4 = Array<T, 8>{
    static_cast<T>(0.032223100604042702),
    static_cast<T>(-0.012603967262037833),
    static_cast<T>(-0.099219543576847216),
    static_cast<T>(0.29785779560527736),
    static_cast<T>(0.80373875180591614),
    static_cast<T>(0.49761866763201545),
    static_cast<T>(-0.02963552764599851),
    static_cast<T>(-0.075765714789273325),
};

template<typename T>
inline constexpr auto sym5 = Array<T, 10>{
    static_cast<T>(0.019538882735286728),
    static_cast<T>(-0.021101834024758855),
    static_cast<T>(-0.17532808990845047),
    static_cast<T>(0.016602105764522319),
    static_cast<T>(0.63397896345821192),
    static_cast<T>(0.72340769040242059),
    static_cast<T>(0.1993975339773936),
    static_cast<T>(-0.039134249302383094),
    static_cast<T>(0.029519490925774643),
    static_cast<T>(0.027333068345077982),
};

template<typename T>
inline constexpr auto sym6 = Array<T, 12>{
    static_cast<T>(-0.007800708325034148),
    static_cast<T>(0.0017677118642428036),
    static_cast<T>(0.044724901770665779),
    static_cast<T>(-0.021060292512300564),
    static_cast<T>(-0.072637522786462516),
    static_cast<T>(0.3379294217276218),
    static_cast<T>(0.787641141030194),
    static_cast<T>(0.49105594192674662),
    static_cast<T>(-0.048311742585632998),
    static_cast<T>(-0.11799011114819057),
    static_cast<T>(0.0034907120842174702),
    static_cast<T>(0.015404109327027373),
};

template<typename T>
inline constexpr auto sym7 = Array<T, 14>{
    static_cast<T>(0.010268176708511255),
    static_cast<T>(0.0040102448715336634),
    static_cast<T>(-0.10780823770381774),
    static_cast<T>(-0.14004724044296152),
    static_cast<T>(0.28862963175151463),
    static_cast<T>(0.76776431700316405),
    static_cast<T>(0.5361019170917628),
    static_cast<T>(0.017441255086855827),
    static_cast<T>(-0.049552834937127255),
    static_cast<T>(0.067892693501372697),
    static_cast<T>(0.03051551316596357),
    static_cast<T>(-0.01263630340325193),
    static_cast<T>(-0.0010473848886829163),
    static_cast<T>(0.0026818145682578781),
};

template<typename T>
inline constexpr auto sym8 = Array<T, 16>{
    static_cast<T>(0.0018899503327594609),
    static_cast<T>(-0.0003029205147213668),
    static_cast<T>(-0.014952258337048231),
    static_cast<T>(0.0038087520138906151),
    static_cast<T>(0.049137179673607506),
    static_cast<T>(-0.027219029917056003),
    static_cast<T>(-0.051945838107709037),
    static_cast<T>(0.3644418948353314),
    static_cast<T>(0.77718575170052351),
    static_cast<T>(0.48135965125837221),
    static_cast<T>(-0.061273359067658524),
    static_cast<T>(-0.14329423835080971),
    static_cast<T>(0.0076074873249176054),
    static_cast<T>(0.031695087811492981),
    static_cast<T>(-0.00054213233179114812),
    static_cast<T>(-0.0033824159510061256),
};

template<typename T>
inline constexpr auto sym9 = Array<T, 18>{
    static_cast<T>(0.0010694900329086053),
    static_cast<T>(-0.00047315449868008311),
    static_cast<T>(-0.010264064027633142),
    static_cast<T>(0.0088592674934004842),
    static_cast<T>(0.06207778930288603),
    static_cast<T>(-0.018233770779395985),
    static_cast<T>(-0.19155083129728512),
    static_cast<T>(0.035272488035271894),
    static_cast<T>(0.61733844914093583),
    static_cast<T>(0.717897082764412),
    static_cast<T>(0.238760914607303),
    static_cast<T>(-0.054568958430834071),
    static_cast<T>(0.00058346274612580684),
    static_cast<T>(0.03022487885827568),
    static_cast<T>(-0.01152821020767923),
    static_cast<T>(-0.013271967781817119),
    static_cast<T>(0.00061978088898558676),
    static_cast<T>(0.0014009155259146807),
};

template<typename T>
inline constexpr auto sym10 = Array<T, 20>{
    static_cast<T>(-0.00045932942100465878), static_cast<T>(5.7036083618494284e-005),
    static_cast<T>(0.0045931735853118284),   static_cast<T>(-0.00080435893201654491),
    static_cast<T>(-0.02035493981231129),    static_cast<T>(0.0057649120335819086),
    static_cast<T>(0.049994972077376687),    static_cast<T>(-0.0319900568824278),
    static_cast<T>(-0.035536740473817552),   static_cast<T>(0.38382676106708546),
    static_cast<T>(0.7695100370211071),      static_cast<T>(0.47169066693843925),
    static_cast<T>(-0.070880535783243853),   static_cast<T>(-0.15949427888491757),
    static_cast<T>(0.011609893903711381),    static_cast<T>(0.045927239231092203),
    static_cast<T>(-0.0014653825813050513),  static_cast<T>(-0.0086412992770224222),
    static_cast<T>(9.5632670722894754e-005), static_cast<T>(0.00077015980911449011),
};

template<typename T>
inline constexpr auto sym11 = Array<T, 22>{
    static_cast<T>(0.00048926361026192387),   static_cast<T>(0.00011053509764272153),
    static_cast<T>(-0.0063896036664548919),   static_cast<T>(-0.0020034719001093887),
    static_cast<T>(0.043000190681552281),     static_cast<T>(0.035266759564466552),
    static_cast<T>(-0.14460234370531561),     static_cast<T>(-0.2046547944958006),
    static_cast<T>(0.23768990904924897),      static_cast<T>(0.73034354908839572),
    static_cast<T>(0.57202297801008706),      static_cast<T>(0.097198394458909473),
    static_cast<T>(-0.022832651022562687),    static_cast<T>(0.069976799610734136),
    static_cast<T>(0.0370374159788594),       static_cast<T>(-0.024080841595864003),
    static_cast<T>(-0.0098579348287897942),   static_cast<T>(0.0065124956747714497),
    static_cast<T>(0.00058835273539699145),   static_cast<T>(-0.0017343662672978692),
    static_cast<T>(-3.8795655736158566e-005), static_cast<T>(0.00017172195069934854),
};

template<typename T>
inline constexpr auto sym12 = Array<T, 24>{
    static_cast<T>(-0.00017906658697508691),  static_cast<T>(-1.8158078862617515e-005),
    static_cast<T>(0.0023502976141834648),    static_cast<T>(0.00030764779631059454),
    static_cast<T>(-0.014589836449234145),    static_cast<T>(-0.0026043910313322326),
    static_cast<T>(0.057804179445505657),     static_cast<T>(0.01530174062247884),
    static_cast<T>(-0.17037069723886492),     static_cast<T>(-0.07833262231634322),
    static_cast<T>(0.46274103121927235),      static_cast<T>(0.76347909778365719),
    static_cast<T>(0.39888597239022),         static_cast<T>(-0.022162306170337816),
    static_cast<T>(-0.035848830736954392),    static_cast<T>(0.049179318299660837),
    static_cast<T>(0.0075537806116804775),    static_cast<T>(-0.024220722675013445),
    static_cast<T>(-0.0014089092443297553),   static_cast<T>(0.007414965517654251),
    static_cast<T>(0.00018021409008538188),   static_cast<T>(-0.0013497557555715387),
    static_cast<T>(-1.1353928041541452e-005), static_cast<T>(0.00011196719424656033),
};

template<typename T>
inline constexpr auto sym13 = Array<T, 26>{
    static_cast<T>(7.0429866906944016e-005),  static_cast<T>(3.6905373423196241e-005),
    static_cast<T>(-0.0007213643851362283),   static_cast<T>(0.00041326119884196064),
    static_cast<T>(0.0056748537601224395),    static_cast<T>(-0.0014924472742598532),
    static_cast<T>(-0.020749686325515677),    static_cast<T>(0.017618296880653084),
    static_cast<T>(0.092926030899137119),     static_cast<T>(0.0088197576704205465),
    static_cast<T>(-0.14049009311363403),     static_cast<T>(0.11023022302137217),
    static_cast<T>(0.64456438390118564),      static_cast<T>(0.69573915056149638),
    static_cast<T>(0.19770481877117801),      static_cast<T>(-0.12436246075153011),
    static_cast<T>(-0.059750627717943698),    static_cast<T>(0.013862497435849205),
    static_cast<T>(-0.017211642726299048),    static_cast<T>(-0.02021676813338983),
    static_cast<T>(0.0052963597387250252),    static_cast<T>(0.0075262253899680996),
    static_cast<T>(-0.00017094285853022211),  static_cast<T>(-0.0011360634389281183),
    static_cast<T>(-3.5738623648689009e-005), static_cast<T>(6.8203252630753188e-005),
};

template<typename T>
inline constexpr auto sym14 = Array<T, 28>{
    static_cast<float>(4.4618977991475265e-005),
    static_cast<float>(1.9329016965523917e-005),
    static_cast<float>(-0.00060576018246643346),
    static_cast<float>(-7.3214213567023991e-005),
    static_cast<float>(0.0045326774719456481),
    static_cast<float>(0.0010131419871842082),
    static_cast<float>(-0.019439314263626713),
    static_cast<float>(-0.0023650488367403851),
    static_cast<float>(0.069827616361807551),
    static_cast<float>(0.025898587531046669),
    static_cast<float>(-0.15999741114652205),
    static_cast<float>(-0.058111823317717831),
    static_cast<float>(0.47533576263420663),
    static_cast<float>(0.75997624196109093),
    static_cast<float>(0.39320152196208885),
    static_cast<float>(-0.035318112114979733),
    static_cast<float>(-0.057634498351326995),
    static_cast<float>(0.037433088362853452),
    static_cast<float>(0.0042805204990193782),
    static_cast<float>(-0.029196217764038187),
    static_cast<float>(-0.0027537747912240711),
    static_cast<float>(0.010037693717672269),
    static_cast<float>(0.00036647657366011829),
    static_cast<float>(-0.002579441725933078),
    static_cast<float>(-6.2865424814776362e-005),
    static_cast<float>(0.00039843567297594335),
    static_cast<float>(1.1210865808890361e-005),
    static_cast<float>(-2.5879090265397886e-005),
};

template<typename T>
inline constexpr auto sym15 = Array<T, 30>{
    static_cast<T>(2.8660708525318081e-005),  static_cast<T>(2.1717890150778919e-005),
    static_cast<T>(-0.00040216853760293483),  static_cast<T>(-0.00010815440168545525),
    static_cast<T>(0.003481028737064895),     static_cast<T>(0.0015261382781819983),
    static_cast<T>(-0.017171252781638731),    static_cast<T>(-0.0087447888864779517),
    static_cast<T>(0.067969829044879179),     static_cast<T>(0.068393310060480245),
    static_cast<T>(-0.13405629845625389),     static_cast<T>(-0.1966263587662373),
    static_cast<T>(0.2439627054321663),       static_cast<T>(0.72184302963618119),
    static_cast<T>(0.57864041521503451),      static_cast<T>(0.11153369514261872),
    static_cast<T>(-0.04108266663538248),     static_cast<T>(0.040735479696810677),
    static_cast<T>(0.021937642719753955),     static_cast<T>(-0.038876716876833493),
    static_cast<T>(-0.019405011430934468),    static_cast<T>(0.010079977087905669),
    static_cast<T>(0.003423450736351241),     static_cast<T>(-0.0035901654473726417),
    static_cast<T>(-0.00026731644647180568),  static_cast<T>(0.0010705672194623959),
    static_cast<T>(5.5122547855586653e-005),  static_cast<T>(-0.00016066186637495343),
    static_cast<T>(-7.3596667989194696e-006), static_cast<T>(9.7124197379633478e-006),
};

template<typename T>
inline constexpr auto sym16 = Array<T, 32>{
    static_cast<T>(-1.0797982104319795e-005), static_cast<T>(-5.3964831793152419e-006),
    static_cast<T>(0.00016545679579108483),   static_cast<T>(3.656592483348223e-005),
    static_cast<T>(-0.0013387206066921965),   static_cast<T>(-0.00022211647621176323),
    static_cast<T>(0.0069377611308027096),    static_cast<T>(0.001359844742484172),
    static_cast<T>(-0.024952758046290123),    static_cast<T>(-0.0035102750683740089),
    static_cast<T>(0.078037852903419913),     static_cast<T>(0.03072113906330156),
    static_cast<T>(-0.15959219218520598),     static_cast<T>(-0.054040601387606135),
    static_cast<T>(0.47534280601152273),      static_cast<T>(0.75652498787569711),
    static_cast<T>(0.39712293362064416),      static_cast<T>(-0.034574228416972504),
    static_cast<T>(-0.066983049070217779),    static_cast<T>(0.032333091610663785),
    static_cast<T>(0.0048692744049046071),    static_cast<T>(-0.031051202843553064),
    static_cast<T>(-0.0031265171722710075),   static_cast<T>(0.012666731659857348),
    static_cast<T>(0.00071821197883178923),   static_cast<T>(-0.0038809122526038786),
    static_cast<T>(-0.0001084456223089688),   static_cast<T>(0.00085235471080470952),
    static_cast<T>(2.8078582128442894e-005),  static_cast<T>(-0.00010943147929529757),
    static_cast<T>(-3.1135564076219692e-006), static_cast<T>(6.2300067012207606e-006),

};

template<typename T>
inline constexpr auto sym17 = Array<T, 34>{
    static_cast<T>(3.7912531943321266e-006),  static_cast<T>(-2.4527163425832999e-006),
    static_cast<T>(-7.6071244056051285e-005), static_cast<T>(2.5207933140828779e-005),
    static_cast<T>(0.0007198270642148971),    static_cast<T>(5.8400428694052584e-005),
    static_cast<T>(-0.0039323252797979023),   static_cast<T>(-0.0019054076898526659),
    static_cast<T>(0.012396988366648726),     static_cast<T>(0.0099529825235095976),
    static_cast<T>(-0.01803889724191924),     static_cast<T>(-0.0072616347509287674),
    static_cast<T>(0.016158808725919346),     static_cast<T>(-0.086070874720733381),
    static_cast<T>(-0.15507600534974825),     static_cast<T>(0.18053958458111286),
    static_cast<T>(0.68148899534492502),      static_cast<T>(0.65071662920454565),
    static_cast<T>(0.14239835041467819),      static_cast<T>(-0.11856693261143636),
    static_cast<T>(0.0172711782105185),       static_cast<T>(0.10475461484223211),
    static_cast<T>(0.017903952214341119),     static_cast<T>(-0.033291383492359328),
    static_cast<T>(-0.0048192128031761478),   static_cast<T>(0.010482366933031529),
    static_cast<T>(0.0008567700701915741),    static_cast<T>(-0.0027416759756816018),
    static_cast<T>(-0.00013864230268045499),  static_cast<T>(0.0004759963802638669),
    static_cast<T>(-1.3506383399901165e-005), static_cast<T>(-6.2937025975541919e-005),
    static_cast<T>(2.7801266938414138e-006),  static_cast<T>(4.297343327345983e-006),
};

template<typename T>
inline constexpr auto sym18 = Array<T, 36>{
    static_cast<T>(-1.5131530692371587e-006), static_cast<T>(7.8472980558317646e-007),
    static_cast<T>(2.9557437620930811e-005),  static_cast<T>(-9.858816030140058e-006),
    static_cast<T>(-0.00026583011024241041),  static_cast<T>(4.7416145183736671e-005),
    static_cast<T>(0.0014280863270832796),    static_cast<T>(-0.00018877623940755607),
    static_cast<T>(-0.0052397896830266083),   static_cast<T>(0.0010877847895956929),
    static_cast<T>(0.015012356344250213),     static_cast<T>(-0.0032607442000749834),
    static_cast<T>(-0.031712684731814537),    static_cast<T>(0.0062779445543116943),
    static_cast<T>(0.028529597039037808),     static_cast<T>(-0.073799207290607169),
    static_cast<T>(-0.032480573290138676),    static_cast<T>(0.40148386057061813),
    static_cast<T>(0.75362914010179283),      static_cast<T>(0.47396905989393956),
    static_cast<T>(-0.052029158983952786),    static_cast<T>(-0.15993814866932407),
    static_cast<T>(0.033995667103947358),     static_cast<T>(0.084219929970386548),
    static_cast<T>(-0.0050770851607570529),   static_cast<T>(-0.030325091089369604),
    static_cast<T>(0.0016429863972782159),    static_cast<T>(0.0095021643909623654),
    static_cast<T>(-0.00041152110923597756),  static_cast<T>(-0.0023138718145060992),
    static_cast<T>(7.0212734590362685e-005),  static_cast<T>(0.00039616840638254753),
    static_cast<T>(-1.4020992577726755e-005), static_cast<T>(-4.5246757874949856e-005),
    static_cast<T>(1.354915761832114e-006),   static_cast<T>(2.6126125564836423e-006),
};

template<typename T>
inline constexpr auto sym19 = Array<T, 38>{
    static_cast<T>(1.7509367995348687e-006),  static_cast<T>(2.0623170632395688e-006),
    static_cast<T>(-2.8151138661550245e-005), static_cast<T>(-1.6821387029373716e-005),
    static_cast<T>(0.00027621877685734072),   static_cast<T>(0.00012930767650701415),
    static_cast<T>(-0.0017049602611649971),   static_cast<T>(-0.00061792232779831076),
    static_cast<T>(0.0082622369555282547),    static_cast<T>(0.0043193518748949689),
    static_cast<T>(-0.027709896931311252),    static_cast<T>(-0.016908234861345205),
    static_cast<T>(0.084072676279245043),     static_cast<T>(0.093630843415897141),
    static_cast<T>(-0.11624173010739675),     static_cast<T>(-0.17659686625203097),
    static_cast<T>(0.25826616923728363),      static_cast<T>(0.71955552571639425),
    static_cast<T>(0.57814494533860505),      static_cast<T>(0.10902582508127781),
    static_cast<T>(-0.067525058040294086),    static_cast<T>(0.0089545911730436242),
    static_cast<T>(0.0070155738571741596),    static_cast<T>(-0.046635983534938946),
    static_cast<T>(-0.022651993378245951),    static_cast<T>(0.015797439295674631),
    static_cast<T>(0.0079684383206133063),    static_cast<T>(-0.005122205002583014),
    static_cast<T>(-0.0011607032572062486),   static_cast<T>(0.0021214250281823303),
    static_cast<T>(0.00015915804768084938),   static_cast<T>(-0.00063576451500433403),
    static_cast<T>(-4.6120396002105868e-005), static_cast<T>(0.0001155392333357879),
    static_cast<T>(8.8733121737292863e-006),  static_cast<T>(-1.1880518269823984e-005),
    static_cast<T>(-6.4636513033459633e-007), static_cast<T>(5.4877327682158382e-007),
};

template<typename T>
inline constexpr auto sym20 = Array<T, 40>{
    static_cast<T>(-6.3291290447763946e-007), static_cast<T>(-3.2567026420174407e-007),
    static_cast<T>(1.22872527779612e-005),    static_cast<T>(4.5254222091516362e-006),
    static_cast<T>(-0.00011739133516291466),  static_cast<T>(-2.6615550335516086e-005),
    static_cast<T>(0.00074761085978205719),   static_cast<T>(0.00012544091723067259),
    static_cast<T>(-0.0034716478028440734),   static_cast<T>(-0.0006111263857992088),
    static_cast<T>(0.012157040948785737),     static_cast<T>(0.0019385970672402002),
    static_cast<T>(-0.035373336756604236),    static_cast<T>(-0.0068437019650692274),
    static_cast<T>(0.088919668028199561),     static_cast<T>(0.036250951653933078),
    static_cast<T>(-0.16057829841525254),     static_cast<T>(-0.051088342921067398),
    static_cast<T>(0.47199147510148703),      static_cast<T>(0.75116272842273002),
    static_cast<T>(0.40583144434845059),      static_cast<T>(-0.029819368880333728),
    static_cast<T>(-0.078994344928398158),    static_cast<T>(0.025579349509413946),
    static_cast<T>(0.0081232283560096815),    static_cast<T>(-0.031629437144957966),
    static_cast<T>(-0.0033138573836233591),   static_cast<T>(0.017004049023390339),
    static_cast<T>(0.0014230873594621453),    static_cast<T>(-0.0066065857990888609),
    static_cast<T>(-0.0003052628317957281),   static_cast<T>(0.0020889947081901982),
    static_cast<T>(7.2159911880740349e-005),  static_cast<T>(-0.00049473109156726548),
    static_cast<T>(-1.928412300645204e-005),  static_cast<T>(7.992967835772481e-005),
    static_cast<T>(3.0256660627369661e-006),  static_cast<T>(-7.919361411976999e-006),
    static_cast<T>(-1.9015675890554106e-007), static_cast<T>(3.695537474835221e-007),
};
}  // namespace mc
