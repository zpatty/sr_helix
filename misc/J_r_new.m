function Jr = J_r_new(L0,d,dL1,dL2,dL3,dc,dx1,dx2,dx3,dy1,dy2,dy3,th)
%J_r_new
%    Jr = J_r_new(L0,D,dL1,dL2,dL3,DC,DX1,DX2,DX3,DY1,DY2,DY3,TH)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    23-Oct-2024 09:36:41

t2 = cos(th);
t3 = sin(th);
t4 = L0+dL1;
t5 = L0+dL2;
t6 = L0+dL3;
t7 = dx1.*2.0;
t8 = dx2.*2.0;
t9 = dy1.*2.0;
t10 = dy2.*2.0;
t11 = dy3.*2.0;
t12 = dx1.^2;
t13 = dx1.^3;
t14 = dx2.^2;
t15 = dx2.^3;
t16 = dx3.^2;
t17 = dx3.^3;
t18 = dy1.^2;
t19 = dy1.^3;
t20 = dy2.^2;
t21 = dy2.^3;
t22 = dy3.^2;
t23 = dy3.^3;
t28 = 1.0./d;
t24 = dx1.*t2;
t25 = dy1.*t2;
t26 = dx1.*t3;
t27 = dy1.*t3;
t30 = t12+t18;
t31 = t14+t20;
t32 = t16+t22;
t29 = -t27;
t33 = t25+t26;
t34 = 1.0./t30;
t36 = 1.0./t31;
t38 = 1.0./t32;
t40 = sqrt(t30);
t41 = sqrt(t31);
t42 = sqrt(t32);
t35 = t34.^2;
t37 = t36.^2;
t39 = t38.^2;
t43 = t24+t29;
t44 = 1.0./t40;
t46 = 1.0./t41;
t48 = 1.0./t42;
t50 = t28.*t40;
t51 = t28.*t41;
t52 = t28.*t42;
t45 = t44.^3;
t47 = t46.^3;
t49 = t48.^3;
t53 = cos(t50);
t54 = cos(t51);
t55 = cos(t52);
t56 = sin(t50);
t57 = sin(t51);
t58 = sin(t52);
t59 = t53-1.0;
t60 = t54-1.0;
t61 = t55-1.0;
t62 = t7.*t53;
t63 = t8.*t54;
t64 = dx3.*t55.*2.0;
t65 = t9.*t53;
t66 = t10.*t54;
t67 = t12.*t53;
t68 = t14.*t54;
t69 = t16.*t55;
t70 = t18.*t53;
t71 = t20.*t54;
t87 = dx1.*t56.*t57;
t88 = dy1.*t56.*t57;
t89 = t13.*t28.*t44.*t56;
t90 = t15.*t28.*t46.*t57;
t91 = t17.*t28.*t48.*t58;
t92 = t19.*t28.*t44.*t56;
t93 = t21.*t28.*t46.*t57;
t96 = dy1.*t4.*t24.*t45.*t56;
t97 = dy1.*t4.*t26.*t45.*t56;
t102 = dx1.*dx2.*dy1.*t53.*t57;
t103 = dx1.*dy1.*dy2.*t53.*t57;
t104 = dx1.*t18.*t28.*t44.*t56;
t105 = dy1.*t12.*t28.*t44.*t56;
t106 = dx2.*t20.*t28.*t46.*t57;
t107 = dy2.*t14.*t28.*t46.*t57;
t108 = dy3.*t16.*t28.*t48.*t58;
t129 = t18.*t24.*t28.*t45.*t56;
t130 = t12.*t25.*t28.*t45.*t56;
t131 = t18.*t26.*t28.*t45.*t56;
t132 = t12.*t27.*t28.*t45.*t56;
t153 = t46.*t53.*t57;
t154 = d.*dx1.*t41.*t53.*t54;
t155 = d.*dy1.*t41.*t53.*t54;
t156 = d.*dx2.*t40.*t56.*t57;
t157 = d.*dy2.*t40.*t56.*t57;
t160 = dx2.*dy2.*t47.*t53.*t57;
t161 = t2.*t44.*t54.*t56;
t162 = t3.*t44.*t54.*t56;
t165 = t14.*t47.*t53.*t57;
t166 = t20.*t47.*t53.*t57;
t169 = dx2.*t50.*t53.*t57;
t170 = dy2.*t50.*t53.*t57;
t171 = dx2.*dy2.*t28.*t36.*t53.*t54;
t177 = dx1.*dx2.*dy2.*t28.*t46.*t54.*t56;
t178 = dx2.*dy1.*dy2.*t28.*t46.*t54.*t56;
t181 = t33.*t44.*t54.*t56;
t186 = dx1.*t33.*t45.*t54.*t56;
t187 = dy1.*t33.*t45.*t54.*t56;
t188 = t40.*t41.*t53.*t54;
t189 = t43.*t44.*t54.*t56;
t190 = dx1.*t40.*t41.*t54.*t56;
t191 = dy1.*t40.*t41.*t54.*t56;
t192 = dx1.*t28.*t33.*t34.*t53.*t54;
t193 = dy1.*t28.*t33.*t34.*t53.*t54;
t196 = dx1.*t43.*t45.*t54.*t56;
t197 = dy1.*t43.*t45.*t54.*t56;
t199 = dx2.*t40.*t46.*t53.*t54;
t200 = dy2.*t40.*t46.*t53.*t54;
t206 = dx1.*t28.*t34.*t43.*t53.*t54;
t208 = dy1.*t28.*t34.*t43.*t53.*t54;
t218 = dx2.*t2.*t44.*t46.*t56.*t57;
t219 = dy2.*t2.*t44.*t46.*t56.*t57;
t220 = dx2.*t3.*t44.*t46.*t56.*t57;
t221 = dy2.*t3.*t44.*t46.*t56.*t57;
t229 = d.*t3.*t5.*t44.*t46.*t56.*t57;
t238 = d.*t2.*t5.*t44.*t46.*t56.*t57;
t255 = t33.*t44.*t46.*t56.*t57;
t277 = dx2.*dy2.*t33.*t44.*t47.*t56.*t57;
t279 = t43.*t44.*t46.*t56.*t57;
t286 = t14.*t33.*t44.*t47.*t56.*t57;
t287 = t20.*t33.*t44.*t47.*t56.*t57;
t300 = dx2.*dy2.*t43.*t44.*t47.*t56.*t57;
t310 = t14.*t43.*t44.*t47.*t56.*t57;
t312 = t20.*t43.*t44.*t47.*t56.*t57;
t72 = t12+t70;
t73 = t18+t67;
t74 = t14+t71;
t75 = t20+t68;
t76 = t22+t69;
t77 = t24.*t34.*t59;
t78 = t25.*t34.*t59;
t79 = t26.*t34.*t59;
t80 = t27.*t34.*t59;
t86 = t29.*t34.*t59;
t94 = d.*t2.*t4.*t34.*t59;
t95 = d.*t3.*t4.*t34.*t59;
t100 = dx2.*t87;
t101 = dy2.*t88;
t113 = t2.*t7.*t18.*t35.*t59;
t114 = dy1.*t7.*t24.*t35.*t59;
t115 = t3.*t7.*t18.*t35.*t59;
t116 = dy1.*t7.*t26.*t35.*t59;
t117 = -t89;
t118 = -t90;
t119 = -t91;
t120 = -t92;
t121 = -t93;
t123 = dx2.*t57.*t67;
t124 = dy2.*t57.*t70;
t133 = -t104;
t134 = -t105;
t135 = -t106;
t136 = -t107;
t137 = -t108;
t138 = d.*t4.*t7.*t25.*t35.*t59;
t143 = d.*t4.*t7.*t27.*t35.*t59;
t158 = dx2.*t153;
t159 = dy2.*t153;
t163 = -t154;
t164 = -t155;
t167 = -t161;
t168 = -t162;
t172 = -t165;
t173 = -t166;
t174 = t28.*t36.*t53.*t68;
t175 = t28.*t36.*t53.*t71;
t176 = -t171;
t179 = dx1.*t28.*t46.*t56.*t68;
t180 = dy1.*t28.*t46.*t56.*t71;
t198 = -t188;
t207 = -t192;
t209 = -t193;
t210 = -t196;
t211 = -t199;
t212 = -t200;
t213 = dx1.*dx2.*t36.*t44.*t56.*t60;
t214 = dx1.*dy2.*t36.*t44.*t56.*t60;
t215 = dx2.*dy1.*t36.*t44.*t56.*t60;
t216 = dx2.*dy2.*t36.*t44.*t56.*t60;
t217 = dy1.*dy2.*t36.*t44.*t56.*t60;
t222 = -t208;
t225 = dx1.*dx2.*dy1.*dy2.*t36.*t45.*t56.*t60;
t231 = dy2.*t28.*t44.*t46.*t87;
t232 = dx2.*t28.*t44.*t46.*t88;
t234 = dx2.*dy2.*t12.*t36.*t45.*t56.*t60;
t235 = dx2.*dy2.*t18.*t36.*t45.*t56.*t60;
t239 = dx1.*dx2.*dy1.*dy2.*t28.*t34.*t36.*t53.*t60;
t243 = dy2.*t7.*t14.*t37.*t44.*t56.*t60;
t244 = dy1.*t8.*t20.*t37.*t44.*t56.*t60;
t249 = dx2.*dy2.*t28.*t34.*t36.*t60.*t67;
t250 = dx2.*dy2.*t28.*t34.*t36.*t60.*t70;
t251 = dx1.*dx2.*t20.*t37.*t44.*t56.*t60.*-2.0;
t252 = dy1.*dy2.*t14.*t37.*t44.*t56.*t60.*-2.0;
t263 = dy2.*t14.*t28.*t44.*t47.*t87;
t264 = dx2.*t20.*t28.*t44.*t47.*t88;
t272 = dx2.*t255;
t273 = dy2.*t255;
t275 = dy2.*t33.*t45.*t46.*t87;
t276 = dx2.*t33.*t45.*t46.*t88;
t290 = dx2.*t279;
t292 = dy2.*t279;
t296 = dy2.*t43.*t45.*t46.*t87;
t297 = dx2.*t43.*t45.*t46.*t88;
t301 = -t277;
t308 = dx2.*dy2.*t28.*t36.*t181;
t311 = -t286;
t313 = -t287;
t320 = t28.*t33.*t36.*t44.*t56.*t68;
t321 = t28.*t33.*t36.*t44.*t56.*t71;
t324 = -t300;
t330 = dx2.*dy2.*t28.*t36.*t189;
t332 = -t312;
t339 = t28.*t36.*t43.*t44.*t56.*t68;
t340 = t28.*t36.*t43.*t44.*t56.*t71;
t81 = dy1.*t77;
t82 = dy1.*t79;
t83 = -t77;
t84 = -t78;
t85 = -t79;
t109 = t2.*t34.*t72;
t110 = t2.*t34.*t73;
t111 = t3.*t34.*t72;
t112 = t3.*t34.*t73;
t122 = -t94;
t126 = t3.*t7.*t35.*t73;
t128 = t3.*t9.*t35.*t73;
t139 = t2.*t7.*t35.*t72;
t141 = t2.*t9.*t35.*t72;
t144 = t26.*t35.*t72.*-2.0;
t145 = t27.*t35.*t72.*-2.0;
t146 = t24.*t35.*t73.*-2.0;
t147 = t25.*t35.*t73.*-2.0;
t148 = t7+t133;
t149 = t8+t135;
t150 = t9+t134;
t151 = t10+t136;
t152 = t11+t137;
t201 = t62+t117;
t202 = t63+t118;
t203 = t64+t119;
t204 = t65+t120;
t205 = t66+t121;
t223 = dy2.*t213;
t224 = dy2.*t215;
t226 = -t214;
t227 = -t215;
t230 = t28.*t44.*t46.*t100;
t233 = t28.*t44.*t46.*t101;
t236 = t36.*t44.*t56.*t74;
t237 = t36.*t44.*t56.*t75;
t253 = dx1.*dy1.*t36.*t45.*t56.*t74;
t254 = dx1.*dy1.*t36.*t45.*t56.*t75;
t256 = t12.*t36.*t45.*t56.*t75;
t257 = t18.*t36.*t45.*t56.*t74;
t259 = dy2.*t7.*t37.*t44.*t56.*t75;
t260 = dy1.*t8.*t37.*t44.*t56.*t74;
t262 = t20.*t28.*t44.*t47.*t100;
t265 = t14.*t28.*t44.*t47.*t101;
t268 = dx1.*dy1.*t28.*t34.*t36.*t53.*t74;
t269 = dx1.*dy1.*t28.*t34.*t36.*t53.*t75;
t270 = dx1.*dx2.*t37.*t44.*t56.*t75.*-2.0;
t271 = dy1.*dy2.*t37.*t44.*t56.*t74.*-2.0;
t274 = t33.*t45.*t46.*t100;
t278 = t33.*t45.*t46.*t101;
t280 = t28.*t34.*t36.*t67.*t75;
t281 = t28.*t34.*t36.*t70.*t74;
t291 = -t272;
t293 = -t273;
t294 = t43.*t45.*t46.*t100;
t298 = -t275;
t299 = -t276;
t302 = t43.*t45.*t46.*t101;
t305 = dx1.*t28.*t33.*t34.*t158;
t306 = dx1.*t28.*t33.*t34.*t159;
t307 = dy1.*t28.*t33.*t34.*t158;
t309 = dy1.*t28.*t33.*t34.*t159;
t314 = t28.*t272;
t315 = t28.*t273;
t323 = -t296;
t327 = dx1.*t28.*t34.*t43.*t158;
t328 = dx1.*t28.*t34.*t43.*t159;
t329 = dy1.*t28.*t34.*t43.*t158;
t331 = dy1.*t28.*t34.*t43.*t159;
t333 = t28.*t290;
t335 = t28.*t292;
t344 = -t330;
t447 = t100+t101+t198;
t464 = t103+t123+t156+t163+t190;
t465 = t102+t124+t157+t164+t191;
t480 = t87+t169+t178+t179+t211;
t481 = t88+t170+t177+t180+t212;
t182 = t2.*t34.*t148;
t183 = t2.*t34.*t150;
t184 = t3.*t34.*t148;
t185 = t3.*t34.*t150;
t240 = dx1.*t237;
t241 = dy1.*t236;
t284 = -t262;
t285 = -t265;
t295 = -t274;
t303 = -t278;
t316 = t2.*t34.*t201;
t317 = t2.*t34.*t204;
t318 = t3.*t34.*t201;
t319 = t3.*t34.*t204;
t322 = -t294;
t325 = t82+t109;
t326 = t81+t112;
t343 = -t329;
t345 = -t331;
t348 = dx1.*t36.*t44.*t56.*t151;
t349 = dy1.*t36.*t44.*t56.*t149;
t352 = dx1.*t36.*t44.*t56.*t202;
t353 = dy1.*t36.*t44.*t56.*t205;
t364 = -t46.*t57.*(t82-t110);
t365 = -t46.*t57.*(t81-t111);
t375 = -dx2.*t36.*t60.*(t81-t111);
t378 = -dy2.*t36.*t60.*(t82-t110);
t379 = -dy2.*t36.*t60.*(t81-t111);
t382 = -dx2.*dy2.*t47.*t57.*(t82-t110);
t388 = dy2.*t46.*t57.*(t81-t111);
t391 = -t14.*t47.*t57.*(t82-t110);
t402 = -dx2.*dy2.*t28.*t36.*t54.*(t81-t111);
t403 = dx2.*t36.*t60.*(t82-t110);
t409 = -t28.*t36.*t71.*(t81-t111);
t410 = t14.*t47.*t57.*(t82-t110);
t411 = t20.*t47.*t57.*(t81-t111);
t412 = -t36.*t74.*(t81-t111);
t414 = -t8.*t20.*t37.*t60.*(t82-t110);
t418 = dx2.*dy2.*t28.*t36.*t54.*(t82-t110);
t428 = t36.*t75.*(t82-t110);
t433 = -dx2.*t20.*t28.*t47.*t57.*(t82-t110);
t437 = dx2.*t37.*t74.*(t81-t111).*2.0;
t438 = dx2.*t37.*t75.*(t82-t110).*2.0;
t439 = dy2.*t37.*t74.*(t81-t111).*2.0;
t440 = dy2.*t37.*t75.*(t82-t110).*2.0;
t452 = -t36.*t205.*(t81-t111);
t194 = -t182;
t195 = -t185;
t337 = -t317;
t338 = -t318;
t350 = -t348;
t351 = -t349;
t354 = t46.*t57.*t325;
t355 = t46.*t57.*t326;
t358 = dx2.*t36.*t60.*t325;
t359 = dx2.*t36.*t60.*t326;
t360 = dy2.*t36.*t60.*t325;
t361 = dy2.*t36.*t60.*t326;
t362 = dx2.*dy2.*t47.*t57.*t325;
t363 = dx2.*dy2.*t47.*t57.*t326;
t368 = dx2.*t364;
t370 = t14.*t47.*t57.*t326;
t371 = t20.*t47.*t57.*t325;
t372 = dx2.*dy2.*t28.*t36.*t54.*t325;
t373 = dx2.*dy2.*t28.*t36.*t54.*t326;
t387 = dy2.*t375;
t389 = t28.*t36.*t68.*t326;
t390 = t28.*t36.*t71.*t325;
t395 = t36.*t74.*t325;
t396 = t36.*t75.*t326;
t397 = t8.*t20.*t37.*t60.*t325;
t398 = t8.*t20.*t37.*t60.*t326;
t399 = dx2.*dy2.*t8.*t37.*t60.*t325;
t400 = dx2.*dy2.*t8.*t37.*t60.*t326;
t407 = dy2.*t403;
t420 = t8.*t37.*t74.*t325;
t421 = t8.*t37.*t75.*t326;
t422 = t10.*t37.*t74.*t325;
t423 = t10.*t37.*t75.*t326;
t425 = dx2.*t20.*t28.*t47.*t57.*t326;
t426 = dy2.*t14.*t28.*t47.*t57.*t325;
t441 = t36.*t149.*t325;
t442 = t36.*t151.*t326;
t448 = t36.*t202.*t326;
t449 = t36.*t205.*t325;
t454 = t158+t224+t240;
t455 = t159+t223+t241;
t456 = t84+t114+t130+t144+t184;
t457 = t85+t115+t131+t147+t183;
t460 = t83+t113+t129+t145+t319;
t461 = t86+t116+t132+t146+t316;
t510 = t153+t172+t174+t217+t252+t270+t285+t352;
t511 = t153+t173+t175+t213+t251+t271+t284+t353;
t514 = t335+t365+t382+t409+t411+t418;
t534 = t279+t332+t340+t403+t414+t433+t439+t452;
t356 = dx2.*t355;
t357 = dy2.*t354;
t366 = dy2.*t358;
t367 = dy2.*t359;
t376 = -t358;
t377 = -t359;
t380 = -t360;
t381 = -t361;
t424 = dx2.*t28.*t371;
t427 = dy2.*t28.*t370;
t444 = -t441;
t446 = -t442;
t451 = -t448;
t453 = -t449;
t458 = t83+t113+t128+t129+t195;
t459 = t86+t116+t132+t139+t194;
t462 = t84+t114+t126+t130+t338;
t463 = t85+t115+t131+t141+t337;
t466 = dy2.*t46.*t57.*t456;
t467 = dx2.*t46.*t57.*t457;
t470 = dx2.*dy2.*t36.*t60.*t456;
t471 = dx2.*dy2.*t36.*t60.*t457;
t476 = t36.*t74.*t456;
t477 = t36.*t75.*t457;
t482 = dx2.*t46.*t57.*t461;
t483 = dy2.*t46.*t57.*t460;
t486 = dx2.*dy2.*t36.*t60.*t460;
t487 = dx2.*dy2.*t36.*t60.*t461;
t492 = t36.*t74.*t460;
t493 = t36.*t75.*t461;
t499 = t189+t368+t388;
t502 = t290+t387+t428;
t503 = t292+t407+t412;
t508 = t160+t176+t226+t243+t260+t263+t351;
t509 = t160+t176+t227+t244+t259+t264+t350;
t468 = dx2.*t46.*t57.*t458;
t469 = dy2.*t46.*t57.*t459;
t472 = -t466;
t473 = -t467;
t474 = dx2.*dy2.*t36.*t60.*t458;
t475 = dx2.*dy2.*t36.*t60.*t459;
t478 = t36.*t75.*t458;
t479 = t36.*t74.*t459;
t484 = dy2.*t46.*t57.*t463;
t485 = dx2.*t46.*t57.*t462;
t488 = dx2.*dy2.*t36.*t60.*t462;
t489 = dx2.*dy2.*t36.*t60.*t463;
t490 = -t486;
t491 = -t487;
t494 = t36.*t75.*t462;
t495 = t36.*t74.*t463;
t496 = t181+t356+t357;
t497 = -t492;
t498 = -t493;
t500 = t291+t366+t396;
t501 = t293+t367+t395;
t520 = t301+t308+t381+t400+t420+t427+t444;
t521 = t301+t308+t376+t397+t423+t424+t446;
t532 = t255+t311+t320+t380+t399+t421+t426+t451;
t533 = t255+t313+t321+t377+t398+t422+t425+t453;
t516 = t167+t187+t209+t468+t484;
t517 = t168+t186+t207+t469+t485;
t518 = t161+t206+t210+t472+t482;
t519 = t162+t197+t222+t473+t483;
t522 = t218+t299+t307+t478+t489;
t523 = t219+t303+t309+t474+t495;
t524 = t220+t295+t305+t475+t494;
t525 = t221+t298+t306+t479+t488;
t526 = t218+t322+t327+t470+t498;
t527 = t219+t323+t328+t476+t491;
t528 = t220+t297+t343+t477+t490;
t529 = t221+t302+t345+t471+t497;
et1 = t5.*t391+dc.*(-t38.*t76.*(t279-t310+t339+t379-t438+t36.*t202.*(t82-t110)+dx2.*dy2.*t8.*t37.*t60.*(t81-t111)+dy2.*t14.*t28.*t47.*t57.*(t81-t111))+dx3.*t48.*t58.*(t333+t391+t402+t46.*t57.*(t82-t110)+t28.*t36.*t68.*(t82-t110)+dx2.*dy2.*t47.*t57.*(t81-t111))+dx3.*dy3.*t38.*t61.*(t300+t344+t378-t437+dy2.*t28.*t410+t36.*t149.*(t81-t111)+dx2.*dy2.*t8.*t37.*t60.*(t82-t110)))+dx2.*t5.*t36.*t189+d.*t5.*t36.*t60.*(t82-t110)-d.*t6.*t48.*t58.*(t333+t391+t402+t46.*t57.*(t82-t110)+t28.*t36.*t68.*(t82-t110)+dx2.*dy2.*t47.*t57.*(t81-t111));
et2 = dx2.*dy2.*t5.*t47.*t57.*(t81-t111)+d.*dx3.*t6.*t38.*t61.*(t279-t310+t339+t379-t438+t36.*t202.*(t82-t110)+dx2.*dy2.*t8.*t37.*t60.*(t81-t111)+dy2.*t14.*t28.*t47.*t57.*(t81-t111))-d.*dy3.*t6.*t38.*t61.*(t300+t344+t378-t437+dy2.*t28.*t410+t36.*t149.*(t81-t111)+dx2.*dy2.*t8.*t37.*t60.*(t82-t110))+d.*dx2.*dy2.*t5.*t37.*t60.*(t81-t111).*2.0-d.*dx2.*t5.*t8.*t37.*t60.*(t82-t110)-d.*dx2.*t5.*t43.*t44.*t47.*t56.*t57;
et3 = dy3.*t41.*t69.*t100+dy3.*t41.*t69.*t101+t23.*t41.*t55.*t100+t23.*t41.*t55.*t101+d.*dy1.*t14.*t16.*t56-d.*dy1.*t14.*t22.*t56-d.*dy1.*t14.*t56.*t69+d.*dy1.*t16.*t56.*t71-d.*dy1.*t22.*t56.*t71-d.*dy1.*t56.*t69.*t71-dy3.*t40.*t53.*t68.*t69-dy3.*t40.*t53.*t69.*t71-t23.*t40.*t53.*t55.*t68-t23.*t40.*t53.*t55.*t71+dy1.*t14.*t22.*t42.*t56.*t58+dy1.*t22.*t42.*t56.*t58.*t71-d.*dx1.*dx2.*dy2.*t16.*t56+d.*dx1.*dx2.*dy2.*t22.*t56-d.*dx1.*dx3.*dy3.*t20.*t56.*2.0+d.*dx1.*dx2.*dy2.*t56.*t69-d.*dx1.*dx3.*dy3.*t56.*t68.*2.0+d.*dy1.*t14.*t22.*t55.*t56+d.*dy1.*t22.*t55.*t56.*t71-d.*dy3.*t41.*t42.*t58.*t100-d.*dy3.*t41.*t42.*t58.*t101+d.*dx3.*dy1.*dy2.*dy3.*t8.*t56+d.*dx1.*dx2.*dy2.*t16.*t54.*t56;
et4 = -d.*dx1.*dx2.*dy2.*t22.*t54.*t56-d.*dx1.*dx2.*dy2.*t22.*t55.*t56-d.*dx1.*dx2.*dy2.*t54.*t56.*t69+d.*dx3.*dy3.*t7.*t20.*t55.*t56+d.*dx3.*dy3.*t7.*t55.*t56.*t68-dx1.*dx2.*dy2.*t22.*t42.*t56.*t58+dx1.*dx3.*dy3.*t20.*t42.*t56.*t58+dx1.*dx3.*dy3.*t42.*t56.*t58.*t68+d.*dy2.*t16.*t40.*t41.*t53.*t57-d.*dy2.*t22.*t40.*t41.*t53.*t57-d.*dy2.*t40.*t41.*t53.*t57.*t69+d.*dy3.*t40.*t42.*t53.*t58.*t68+d.*dy3.*t40.*t42.*t53.*t58.*t71-d.*dx2.*dx3.*dy1.*dy2.*dy3.*t54.*t56.*2.0-d.*dx2.*dx3.*dy1.*dy2.*dy3.*t55.*t56.*2.0+d.*dx3.*dy1.*dy2.*dy3.*t55.*t56.*t63-dx2.*dx3.*dy1.*dy2.*dy3.*t42.*t56.*t58+d.*dx1.*dx2.*dy2.*t22.*t54.*t55.*t56-d.*dx2.*dx3.*dy3.*t40.*t41.*t53.*t57.*2.0+dx1.*dx2.*dy2.*t22.*t42.*t54.*t56.*t58+d.*dy2.*t22.*t40.*t41.*t53.*t55.*t57;
et5 = dy2.*t22.*t40.*t41.*t42.*t53.*t57.*t58+dx2.*dx3.*dy1.*dy2.*dy3.*t42.*t54.*t56.*t58+d.*dx3.*dy3.*t8.*t40.*t41.*t53.*t55.*t57+dx2.*dx3.*dy3.*t40.*t41.*t42.*t53.*t57.*t58;
mt1 = [-dc.*(t38.*t76.*t500-dx3.*t48.*t58.*t496+dx3.*dy3.*t38.*t61.*t501)+d.*t4.*t78+d.*t4.*t79-d.*t5.*t255+d.*t5.*t359+d.*t5.*t360-d.*t6.*t48.*t58.*t496+d.*dx3.*t6.*t38.*t61.*t500+d.*dy3.*t6.*t38.*t61.*t501,-dc.*(t38.*t76.*t502+dx3.*t48.*t58.*t499+dx3.*dy3.*t38.*t61.*t503)+d.*t4.*t80+d.*t4.*t83+d.*t5.*t279+d.*t5.*t379+d.*t5.*t403+d.*t6.*t48.*t58.*t499+d.*dx3.*t6.*t38.*t61.*t502+d.*dy3.*t6.*t38.*t61.*t503,0.0];
mt2 = [-t97+t122+t238-dc.*(t38.*t76.*t526+dx3.*t48.*t58.*t518+dx3.*dy3.*t38.*t61.*t527)+d.*t6.*t48.*t58.*t518+dx1.*t5.*t34.*t43.*t153+t2.*t4.*t12.*t45.*t56+d.*t4.*t7.*t24.*t35.*t59-d.*t5.*t43.*t45.*t46.*t87-d.*dx2.*t5.*t36.*t60.*t461+d.*dx3.*t6.*t38.*t61.*t526-d.*dy1.*t4.*t26.*t35.*t59.*2.0+d.*dy2.*t5.*t36.*t60.*t456+d.*dy3.*t6.*t38.*t61.*t527,-t95+t96+t138+t229-dc.*(t38.*t76.*t524-dx3.*t48.*t58.*t517+dx3.*dy3.*t38.*t61.*t525)-d.*t6.*t48.*t58.*t517+dx1.*t5.*t33.*t34.*t153+t3.*t4.*t12.*t45.*t56+d.*t4.*t7.*t26.*t35.*t59-d.*t5.*t33.*t45.*t46.*t87+d.*dx2.*t5.*t36.*t60.*t462+d.*dx3.*t6.*t38.*t61.*t524+d.*dy2.*t5.*t36.*t60.*t459+d.*dy3.*t6.*t38.*t61.*t525];
mt3 = [dc.*(t38.*t76.*(t225+t230-t237-t239+t256-t280)-dx3.*dy3.*t38.*t61.*(t216-t231-t234+t249-t253+t268)-dx1.*dx3.*t45.*t46.*t48.*t58.*t447+dx3.*t28.*t34.*t46.*t48.*t58.*t464)+dx1.*t4.*t34.*t53-t5.*t44.*t46.*t87-d.*dx1.*t4.*t45.*t56+dx2.*t5.*t34.*t36.*t60.*t67-t6.*t34.*t46.*t48.*t58.*t464-d.*dx3.*t6.*t38.*t61.*(t225+t230-t237-t239+t256-t280)+d.*dy3.*t6.*t38.*t61.*(t216-t231-t234+t249-t253+t268)+d.*dx2.*t5.*t36.*t44.*t56.*t60+dx1.*dy1.*dy2.*t5.*t34.*t36.*t53.*t60-d.*dx2.*t5.*t12.*t36.*t45.*t56.*t60+d.*dx1.*t6.*t45.*t46.*t48.*t58.*t447-d.*dx1.*dy1.*dy2.*t5.*t36.*t45.*t56.*t60];
mt4 = [t95+t96+t138-t229+dc.*(t38.*t76.*t528+dx3.*t48.*t58.*t519+dx3.*dy3.*t38.*t61.*t529)-d.*t6.*t48.*t58.*t519+dy1.*t5.*t34.*t43.*t153-t3.*t4.*t18.*t45.*t56-d.*t3.*t4.*t18.*t35.*t59.*2.0-d.*t5.*t43.*t45.*t46.*t88-d.*dx2.*t5.*t36.*t60.*t457-d.*dx3.*t6.*t38.*t61.*t528+d.*dy2.*t5.*t36.*t60.*t460-d.*dy3.*t6.*t38.*t61.*t529,t97+t122+t143+t238-dc.*(t38.*t76.*t522-dx3.*t48.*t58.*t516+dx3.*dy3.*t38.*t61.*t523)-d.*t6.*t48.*t58.*t516+dy1.*t5.*t33.*t34.*t153+t2.*t4.*t18.*t45.*t56+d.*t4.*t9.*t25.*t35.*t59-d.*t5.*t33.*t45.*t46.*t88+d.*dx2.*t5.*t36.*t60.*t458+d.*dx3.*t6.*t38.*t61.*t522+d.*dy2.*t5.*t36.*t60.*t463+d.*dy3.*t6.*t38.*t61.*t523];
mt5 = [-dc.*(t38.*t76.*(t216-t232-t235+t250-t254+t269)-dx3.*dy3.*t38.*t61.*(t225+t233-t236-t239+t257-t281)+dx3.*dy1.*t45.*t46.*t48.*t58.*t447-dx3.*t28.*t34.*t46.*t48.*t58.*t465)+dy1.*t4.*t34.*t53-t5.*t44.*t46.*t88-d.*dy1.*t4.*t45.*t56+dy2.*t5.*t34.*t36.*t60.*t70-t6.*t34.*t46.*t48.*t58.*t465+d.*dx3.*t6.*t38.*t61.*(t216-t232-t235+t250-t254+t269)-d.*dy3.*t6.*t38.*t61.*(t225+t233-t236-t239+t257-t281)+d.*dy2.*t5.*t36.*t44.*t56.*t60+dx1.*dx2.*dy1.*t5.*t34.*t36.*t53.*t60-d.*dy2.*t5.*t18.*t36.*t45.*t56.*t60+d.*dy1.*t6.*t45.*t46.*t48.*t58.*t447-d.*dx1.*dx2.*dy1.*t5.*t36.*t45.*t56.*t60];
mt6 = [-d.*t34.*t43.*t59,-d.*t33.*t34.*t59,d.*t44.*t56,et1+et2,-dc.*(t38.*t76.*t532-dx3.*t48.*t58.*(t314-t355+t362+t370-t372-t389)+dx3.*dy3.*t38.*t61.*t520)+t5.*t362+t5.*t370+dx2.*t5.*t36.*t181-d.*t6.*t48.*t58.*(t314-t355+t362+t370-t372-t389)-d.*t5.*t36.*t60.*t326+d.*dx3.*t6.*t38.*t61.*t532+d.*dy3.*t6.*t38.*t61.*t520+d.*dx2.*t5.*t8.*t37.*t60.*t326+d.*dy2.*t5.*t8.*t37.*t60.*t325-d.*dx2.*t5.*t33.*t44.*t47.*t56.*t57];
mt7 = [-dc.*(t38.*t76.*t510-dx3.*dy3.*t38.*t61.*t508-dx3.*t44.*t46.*t48.*t58.*t480+dx2.*dx3.*t44.*t47.*t48.*t58.*t447)+dx2.*t5.*t36.*t53.*t54-dx2.*t5.*t44.*t47.*t101-t5.*t14.*t44.*t47.*t87-d.*dx2.*t5.*t47.*t53.*t57+d.*dx3.*t6.*t38.*t61.*t510-d.*dy3.*t6.*t38.*t61.*t508+d.*dx1.*t5.*t36.*t44.*t56.*t60-d.*t6.*t44.*t46.*t48.*t58.*t480-d.*dx1.*t5.*t14.*t37.*t44.*t56.*t60.*2.0+d.*dx2.*t6.*t44.*t47.*t48.*t58.*t447-d.*dx2.*dy1.*dy2.*t5.*t37.*t44.*t56.*t60.*2.0];
mt8 = [-dc.*(t38.*t76.*(t324+t330+t375-t440+dx2.*t28.*t411+t36.*t151.*(t82-t110)+t8.*t20.*t37.*t60.*(t81-t111))-dx3.*t48.*t58.*t514+dx3.*dy3.*t38.*t61.*t534)+t5.*t382+t5.*t411+dy2.*t5.*t36.*t189-d.*t6.*t48.*t58.*t514-d.*t5.*t36.*t60.*(t81-t111)+d.*t5.*t20.*t37.*t60.*(t81-t111).*2.0+d.*dx3.*t6.*t38.*t61.*(t324+t330+t375-t440+dx2.*t28.*t411+t36.*t151.*(t82-t110)+t8.*t20.*t37.*t60.*(t81-t111))+d.*dy3.*t6.*t38.*t61.*t534-d.*dy2.*t5.*t8.*t37.*t60.*(t82-t110)-d.*dy2.*t5.*t43.*t44.*t47.*t56.*t57];
mt9 = [-dc.*(t38.*t76.*t521-dx3.*t48.*t58.*(t315-t354+t363+t371-t373-t390)+dx3.*dy3.*t38.*t61.*t533)+t5.*t363+t5.*t371+dy2.*t5.*t36.*t181-d.*t6.*t48.*t58.*(t315-t354+t363+t371-t373-t390)-d.*t5.*t36.*t60.*t325+d.*dx3.*t6.*t38.*t61.*t521+d.*dy3.*t6.*t38.*t61.*t533+d.*dy2.*t5.*t8.*t37.*t60.*t326+d.*dy2.*t5.*t10.*t37.*t60.*t325-d.*dy2.*t5.*t33.*t44.*t47.*t56.*t57];
mt10 = [dc.*(t38.*t76.*t509-dx3.*dy3.*t38.*t61.*t511+dx3.*t44.*t46.*t48.*t58.*t481-dx3.*dy2.*t44.*t47.*t48.*t58.*t447)+dy2.*t5.*t36.*t53.*t54-dy2.*t5.*t44.*t47.*t100-t5.*t20.*t44.*t47.*t88-d.*dx3.*t6.*t38.*t61.*t509-d.*dy2.*t5.*t47.*t53.*t57+d.*dy3.*t6.*t38.*t61.*t511+d.*dy1.*t5.*t36.*t44.*t56.*t60-d.*t6.*t44.*t46.*t48.*t58.*t481-d.*dy1.*t5.*t20.*t37.*t44.*t56.*t60.*2.0+d.*dy2.*t6.*t44.*t47.*t48.*t58.*t447-d.*dx1.*dx2.*dy2.*t5.*t37.*t44.*t56.*t60.*2.0,d.*t279+d.*t379+d.*t403,d.*t255+d.*t377+d.*t380,d.*t153+d.*t213+d.*t217];
mt11 = [-dc.*(t48.*t58.*t499+t38.*t203.*t502-dx3.*t39.*t76.*t502.*2.0+dy3.*t38.*t61.*t503-t16.*t49.*t58.*t499+t28.*t38.*t69.*t499-dy3.*t16.*t39.*t61.*t503.*2.0-dy3.*t16.*t28.*t49.*t58.*t503)+d.*t6.*t38.*t61.*t502+dx3.*t6.*t38.*t55.*t499-t6.*t16.*t49.*t58.*t502-d.*t6.*t16.*t39.*t61.*t502.*2.0-d.*dx3.*t6.*t49.*t58.*t499-dx3.*dy3.*t6.*t49.*t58.*t503-d.*dx3.*dy3.*t6.*t39.*t61.*t503.*2.0];
mt12 = [-dc.*(t48.*t58.*t496-t38.*t203.*t500+dx3.*t39.*t76.*t500.*2.0-dy3.*t38.*t61.*t501-t16.*t49.*t58.*t496+t28.*t38.*t69.*t496+t11.*t16.*t39.*t61.*t501+dy3.*t16.*t28.*t49.*t58.*t501)-d.*t6.*t38.*t61.*t500+dx3.*t6.*t38.*t55.*t496+t6.*t16.*t49.*t58.*t500+d.*t6.*t16.*t39.*t61.*t500.*2.0-d.*dx3.*t6.*t49.*t58.*t496+dx3.*dy3.*t6.*t49.*t58.*t501+d.*dx3.*t6.*t11.*t39.*t61.*t501];
mt13 = [dc.*(-t38.*t203.*t454+dx3.*t39.*t76.*t454.*2.0-dy3.*t38.*t61.*t455+t11.*t16.*t39.*t61.*t455+t44.*t46.*t48.*t58.*t447+dy3.*t16.*t28.*t49.*t58.*t455-t16.*t44.*t46.*t49.*t58.*t447+t28.*t38.*t44.*t46.*t69.*t447)+d.*t6.*t38.*t61.*t454-t6.*t16.*t49.*t58.*t454-d.*t6.*t16.*t39.*t61.*t454.*2.0-dx3.*dy3.*t6.*t49.*t58.*t455-d.*dx3.*dy3.*t6.*t39.*t61.*t455.*2.0-dx3.*t6.*t38.*t44.*t46.*t55.*t447+d.*dx3.*t6.*t44.*t46.*t49.*t58.*t447];
mt14 = [dc.*(-t38.*t152.*t502-dx3.*t38.*t61.*t503+t11.*t39.*t76.*t502+dx3.*dy3.*t49.*t58.*t499+dx3.*t22.*t28.*t49.*t58.*t503+dx3.*dy3.*t11.*t39.*t61.*t503-dx3.*dy3.*t28.*t38.*t55.*t499)+d.*t6.*t38.*t61.*t503+dy3.*t6.*t38.*t55.*t499-t6.*t22.*t49.*t58.*t503-d.*t6.*t22.*t39.*t61.*t503.*2.0-d.*dy3.*t6.*t49.*t58.*t499-dx3.*dy3.*t6.*t49.*t58.*t502-d.*dx3.*dy3.*t6.*t39.*t61.*t502.*2.0];
mt15 = [-dc.*(-t38.*t152.*t500-dx3.*t38.*t61.*t501+t11.*t39.*t76.*t500-dx3.*dy3.*t49.*t58.*t496+dx3.*t22.*t28.*t49.*t58.*t501+dx3.*dy3.*t11.*t39.*t61.*t501+dx3.*dy3.*t28.*t38.*t55.*t496)-d.*t6.*t38.*t61.*t501+dy3.*t6.*t38.*t55.*t496+t6.*t22.*t49.*t58.*t501-d.*dy3.*t6.*t49.*t58.*t496+dx3.*dy3.*t6.*t49.*t58.*t500+d.*dx3.*t6.*t11.*t39.*t61.*t500+d.*dy3.*t6.*t11.*t39.*t61.*t501,-t28.*t36.*t39.*t44.*(et3+et4+et5).*(L0.*d+d.*dL3-dc.*dx3),d.*t48.*t58.*t499+d.*dx3.*t38.*t61.*t502+d.*dy3.*t38.*t61.*t503,d.*t48.*t58.*t496-d.*dx3.*t38.*t61.*t500-d.*dy3.*t38.*t61.*t501,d.*dx3.*t38.*t61.*t454+d.*dy3.*t38.*t61.*t455-d.*t44.*t46.*t48.*t58.*t447];
Jr = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12,mt13,mt14,mt15],3,10);
end
