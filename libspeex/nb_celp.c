/* Copyright (C) 2002-2006 Jean-Marc Valin 
   File: nb_celp.c

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   
   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   
   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   
   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
   
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS&#39;&#39; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
#include &quot;config.h&quot;
#endif

#include &lt;math.h&gt;
#include &quot;nb_celp.h&quot;
#include &quot;lpc.h&quot;
#include &quot;lsp.h&quot;
#include &quot;ltp.h&quot;
#include &quot;quant_lsp.h&quot;
#include &quot;cb_search.h&quot;
#include &quot;filters.h&quot;
#include &quot;stack_alloc.h&quot;
#include &quot;vq.h&quot;
#include &lt;speex/speex_bits.h&gt;
#include &quot;vbr.h&quot;
#include &quot;arch.h&quot;
#include &quot;math_approx.h&quot;
#include &quot;os_support.h&quot;
#include &lt;speex/speex_callbacks.h&gt;

#ifdef VORBIS_PSYCHO
#include &quot;vorbis_psy.h&quot;
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef NULL
#define NULL 0
#endif

#define SUBMODE(x) st-&gt;submodes[st-&gt;submodeID]-&gt;x

/* Default size for the encoder and decoder stack (can be changed at compile time).
   This does not apply when using variable-size arrays or alloca. */
#ifndef NB_ENC_STACK
#define NB_ENC_STACK (8000*sizeof(spx_sig_t))
#endif

#ifndef NB_DEC_STACK
#define NB_DEC_STACK (4000*sizeof(spx_sig_t))
#endif


#ifdef FIXED_POINT
const spx_word32_t ol_gain_table[32]={18900, 25150, 33468, 44536, 59265, 78865, 104946, 139653, 185838, 247297, 329081, 437913, 582736, 775454, 1031906, 1373169, 1827293, 2431601, 3235761, 4305867, 5729870, 7624808, 10146425, 13501971, 17967238, 23909222, 31816294, 42338330, 56340132, 74972501, 99766822, 132760927};
const spx_word16_t exc_gain_quant_scal3_bound[7]={1841, 3883, 6051, 8062, 10444, 13580, 18560};
const spx_word16_t exc_gain_quant_scal3[8]={1002, 2680, 5086, 7016, 9108, 11781, 15380, 21740};
const spx_word16_t exc_gain_quant_scal1_bound[1]={14385};
const spx_word16_t exc_gain_quant_scal1[2]={11546, 17224};

#define LSP_MARGIN 16
#define LSP_DELTA1 6553
#define LSP_DELTA2 1638

#else

const float exc_gain_quant_scal3_bound[7]={0.112338f, 0.236980f, 0.369316f, 0.492054f, 0.637471f, 0.828874f, 1.132784f};
const float exc_gain_quant_scal3[8]={0.061130f, 0.163546f, 0.310413f, 0.428220f, 0.555887f, 0.719055f, 0.938694f, 1.326874f};
const float exc_gain_quant_scal1_bound[1]={0.87798f};
const float exc_gain_quant_scal1[2]={0.70469f, 1.05127f};

#define LSP_MARGIN .002f
#define LSP_DELTA1 .2f
#define LSP_DELTA2 .05f

#endif

#ifdef VORBIS_PSYCHO
#define EXTRA_BUFFER 100
#else
#define EXTRA_BUFFER 0
#endif


#define sqr(x) ((x)*(x))

extern const spx_word16_t lag_window[];
extern const spx_word16_t lpc_window[];

void *nb_encoder_init(const SpeexMode *m)
{
   EncState *st;
   const SpeexNBMode *mode;
   int i;

   mode=(const SpeexNBMode *)m-&gt;mode;
   st = (EncState*)speex_alloc(sizeof(EncState));
   if (!st)
      return NULL;
#if defined(VAR_ARRAYS) || defined (USE_ALLOCA)
   st-&gt;stack = NULL;
#else
   st-&gt;stack = (char*)speex_alloc_scratch(NB_ENC_STACK);
#endif
   
   st-&gt;mode=m;

   st-&gt;frameSize = mode-&gt;frameSize;
   st-&gt;nbSubframes=mode-&gt;frameSize/mode-&gt;subframeSize;
   st-&gt;subframeSize=mode-&gt;subframeSize;
   st-&gt;windowSize = st-&gt;frameSize+st-&gt;subframeSize;
   st-&gt;lpcSize = mode-&gt;lpcSize;
   st-&gt;gamma1=mode-&gt;gamma1;
   st-&gt;gamma2=mode-&gt;gamma2;
   st-&gt;min_pitch=mode-&gt;pitchStart;
   st-&gt;max_pitch=mode-&gt;pitchEnd;
   st-&gt;lpc_floor = mode-&gt;lpc_floor;
  
   st-&gt;submodes=mode-&gt;submodes;
   st-&gt;submodeID=st-&gt;submodeSelect=mode-&gt;defaultSubmode;
   st-&gt;bounded_pitch = 1;

   st-&gt;encode_submode = 1;

#ifdef VORBIS_PSYCHO
   st-&gt;psy = vorbis_psy_init(8000, 256);
   st-&gt;curve = (float*)speex_alloc(128*sizeof(float));
   st-&gt;old_curve = (float*)speex_alloc(128*sizeof(float));
   st-&gt;psy_window = (float*)speex_alloc(256*sizeof(float));
#endif

   st-&gt;cumul_gain = 1024;

   /* Allocating input buffer */
   st-&gt;winBuf = (spx_word16_t*)speex_alloc((st-&gt;windowSize-st-&gt;frameSize)*sizeof(spx_word16_t));
   /* Allocating excitation buffer */
   st-&gt;excBuf = (spx_word16_t*)speex_alloc((mode-&gt;frameSize+mode-&gt;pitchEnd+2)*sizeof(spx_word16_t));
   st-&gt;exc = st-&gt;excBuf + mode-&gt;pitchEnd + 2;
   st-&gt;swBuf = (spx_word16_t*)speex_alloc((mode-&gt;frameSize+mode-&gt;pitchEnd+2)*sizeof(spx_word16_t));
   st-&gt;sw = st-&gt;swBuf + mode-&gt;pitchEnd + 2;

   st-&gt;window= lpc_window;
   
   /* Create the window for autocorrelation (lag-windowing) */
   st-&gt;lagWindow = lag_window;

   st-&gt;old_lsp = (spx_lsp_t*)speex_alloc((st-&gt;lpcSize)*sizeof(spx_lsp_t));
   st-&gt;old_qlsp = (spx_lsp_t*)speex_alloc((st-&gt;lpcSize)*sizeof(spx_lsp_t));
   st-&gt;first = 1;
   for (i=0;i&lt;st-&gt;lpcSize;i++)
      st-&gt;old_lsp[i]= DIV32(MULT16_16(QCONST16(3.1415927f, LSP_SHIFT), i+1), st-&gt;lpcSize+1);

   st-&gt;mem_sp = (spx_mem_t*)speex_alloc((st-&gt;lpcSize)*sizeof(spx_mem_t));
   st-&gt;mem_sw = (spx_mem_t*)speex_alloc((st-&gt;lpcSize)*sizeof(spx_mem_t));
   st-&gt;mem_sw_whole = (spx_mem_t*)speex_alloc((st-&gt;lpcSize)*sizeof(spx_mem_t));
   st-&gt;mem_exc = (spx_mem_t*)speex_alloc((st-&gt;lpcSize)*sizeof(spx_mem_t));
   st-&gt;mem_exc2 = (spx_mem_t*)speex_alloc((st-&gt;lpcSize)*sizeof(spx_mem_t));

   st-&gt;pi_gain = (spx_word32_t*)speex_alloc((st-&gt;nbSubframes)*sizeof(spx_word32_t));
   st-&gt;innov_rms_save = NULL;
   
   st-&gt;pitch = (int*)speex_alloc((st-&gt;nbSubframes)*sizeof(int));

#ifndef DISABLE_VBR
   st-&gt;vbr = (VBRState*)speex_alloc(sizeof(VBRState));
   vbr_init(st-&gt;vbr);
   st-&gt;vbr_quality = 8;
   st-&gt;vbr_enabled = 0;
   st-&gt;vbr_max = 0;
   st-&gt;vad_enabled = 0;
   st-&gt;dtx_enabled = 0;
   st-&gt;dtx_count=0;
   st-&gt;abr_enabled = 0;
   st-&gt;abr_drift = 0;
   st-&gt;abr_drift2 = 0;
#endif /* #ifndef DISABLE_VBR */

   st-&gt;plc_tuning = 2;
   st-&gt;complexity=2;
   st-&gt;sampling_rate=8000;
   st-&gt;isWideband = 0;
   st-&gt;highpass_enabled = 1;
   
#ifdef ENABLE_VALGRIND
   VALGRIND_MAKE_READABLE(st, NB_ENC_STACK);
#endif
   return st;
}

void nb_encoder_destroy(void *state)
{
   EncState *st=(EncState *)state;
   /* Free all allocated memory */
#if !(defined(VAR_ARRAYS) || defined (USE_ALLOCA))
   speex_free_scratch(st-&gt;stack);
#endif

   speex_free (st-&gt;winBuf);
   speex_free (st-&gt;excBuf);
   speex_free (st-&gt;old_qlsp);
   speex_free (st-&gt;swBuf);

   speex_free (st-&gt;old_lsp);
   speex_free (st-&gt;mem_sp);
   speex_free (st-&gt;mem_sw);
   speex_free (st-&gt;mem_sw_whole);
   speex_free (st-&gt;mem_exc);
   speex_free (st-&gt;mem_exc2);
   speex_free (st-&gt;pi_gain);
   speex_free (st-&gt;pitch);

#ifndef DISABLE_VBR
   vbr_destroy(st-&gt;vbr);
   speex_free (st-&gt;vbr);
#endif /* #ifndef DISABLE_VBR */

#ifdef VORBIS_PSYCHO
   vorbis_psy_destroy(st-&gt;psy);
   speex_free (st-&gt;curve);
   speex_free (st-&gt;old_curve);
   speex_free (st-&gt;psy_window);
#endif

   /*Free state memory... should be last*/
   speex_free(st);
}

int nb_encode(void *state, void *vin, SpeexBits *bits)
{
   EncState *st;
   int i, sub, roots;
   int ol_pitch;
   spx_word16_t ol_pitch_coef;
   spx_word32_t ol_gain;
   VARDECL(spx_word16_t *ringing);
   VARDECL(spx_word16_t *target);
   VARDECL(spx_sig_t *innov);
   VARDECL(spx_word32_t *exc32);
   VARDECL(spx_mem_t *mem);
   VARDECL(spx_coef_t *bw_lpc1);
   VARDECL(spx_coef_t *bw_lpc2);
   VARDECL(spx_coef_t *lpc);
   VARDECL(spx_lsp_t *lsp);
   VARDECL(spx_lsp_t *qlsp);
   VARDECL(spx_lsp_t *interp_lsp);
   VARDECL(spx_lsp_t *interp_qlsp);
   VARDECL(spx_coef_t *interp_lpc);
   VARDECL(spx_coef_t *interp_qlpc);
   char *stack;
   VARDECL(spx_word16_t *syn_resp);
   VARDECL(spx_word16_t *real_exc);
   
   spx_word32_t ener=0;
   spx_word16_t fine_gain;
   spx_word16_t *in = (spx_word16_t*)vin;

   st=(EncState *)state;
   stack=st-&gt;stack;

   ALLOC(lpc, st-&gt;lpcSize, spx_coef_t);
   ALLOC(bw_lpc1, st-&gt;lpcSize, spx_coef_t);
   ALLOC(bw_lpc2, st-&gt;lpcSize, spx_coef_t);
   ALLOC(lsp, st-&gt;lpcSize, spx_lsp_t);
   ALLOC(qlsp, st-&gt;lpcSize, spx_lsp_t);
   ALLOC(interp_lsp, st-&gt;lpcSize, spx_lsp_t);
   ALLOC(interp_qlsp, st-&gt;lpcSize, spx_lsp_t);
   ALLOC(interp_lpc, st-&gt;lpcSize, spx_coef_t);
   ALLOC(interp_qlpc, st-&gt;lpcSize, spx_coef_t);

   /* Move signals 1 frame towards the past */
   SPEEX_MOVE(st-&gt;excBuf, st-&gt;excBuf+st-&gt;frameSize, st-&gt;max_pitch+2);
   SPEEX_MOVE(st-&gt;swBuf, st-&gt;swBuf+st-&gt;frameSize, st-&gt;max_pitch+2);

   if (st-&gt;highpass_enabled)
      highpass(in, in, st-&gt;frameSize, (st-&gt;isWideband?HIGHPASS_WIDEBAND:HIGHPASS_NARROWBAND)|HIGHPASS_INPUT, st-&gt;mem_hp);
   
   {
      VARDECL(spx_word16_t *w_sig);
      VARDECL(spx_word16_t *autocorr);
      ALLOC(w_sig, st-&gt;windowSize, spx_word16_t);
      ALLOC(autocorr, st-&gt;lpcSize+1, spx_word16_t);
      /* Window for analysis */
      for (i=0;i&lt;st-&gt;windowSize-st-&gt;frameSize;i++)
         w_sig[i] = EXTRACT16(SHR32(MULT16_16(st-&gt;winBuf[i],st-&gt;window[i]),SIG_SHIFT));
      for (;i&lt;st-&gt;windowSize;i++)
         w_sig[i] = EXTRACT16(SHR32(MULT16_16(in[i-st-&gt;windowSize+st-&gt;frameSize],st-&gt;window[i]),SIG_SHIFT));
      /* Compute auto-correlation */
      _spx_autocorr(w_sig, autocorr, st-&gt;lpcSize+1, st-&gt;windowSize);
      autocorr[0] = ADD16(autocorr[0],MULT16_16_Q15(autocorr[0],st-&gt;lpc_floor)); /* Noise floor in auto-correlation domain */

      /* Lag windowing: equivalent to filtering in the power-spectrum domain */
      for (i=0;i&lt;st-&gt;lpcSize+1;i++)
         autocorr[i] = MULT16_16_Q14(autocorr[i],st-&gt;lagWindow[i]);

      /* Levinson-Durbin */
      _spx_lpc(lpc, autocorr, st-&gt;lpcSize);
      /* LPC to LSPs (x-domain) transform */
      roots=lpc_to_lsp (lpc, st-&gt;lpcSize, lsp, 10, LSP_DELTA1, stack);
      /* Check if we found all the roots */
      if (roots!=st-&gt;lpcSize)
      {
         /*If we can&#39;t find all LSP&#39;s, do some damage control and use previous filter*/
         for (i=0;i&lt;st-&gt;lpcSize;i++)
         {
            lsp[i]=st-&gt;old_lsp[i];
         }
      }
   }




   /* Whole frame analysis (open-loop estimation of pitch and excitation gain) */
   {
      int diff = st-&gt;windowSize-st-&gt;frameSize;
      if (st-&gt;first)
         for (i=0;i&lt;st-&gt;lpcSize;i++)
            interp_lsp[i] = lsp[i];
      else
         lsp_interpolate(st-&gt;old_lsp, lsp, interp_lsp, st-&gt;lpcSize, st-&gt;nbSubframes, st-&gt;nbSubframes&lt;&lt;1);

      lsp_enforce_margin(interp_lsp, st-&gt;lpcSize, LSP_MARGIN);

      /* Compute interpolated LPCs (unquantized) for whole frame*/
      lsp_to_lpc(interp_lsp, interp_lpc, st-&gt;lpcSize,stack);


      /*Open-loop pitch*/
      if (!st-&gt;submodes[st-&gt;submodeID] || (st-&gt;complexity&gt;2 &amp;&amp; SUBMODE(have_subframe_gain)&lt;3) || SUBMODE(forced_pitch_gain) || SUBMODE(lbr_pitch) != -1 
#ifndef DISABLE_VBR
           || st-&gt;vbr_enabled || st-&gt;vad_enabled
#endif
                  )
      {
         int nol_pitch[6];
         spx_word16_t nol_pitch_coef[6];
         
         bw_lpc(st-&gt;gamma1, interp_lpc, bw_lpc1, st-&gt;lpcSize);
         bw_lpc(st-&gt;gamma2, interp_lpc, bw_lpc2, st-&gt;lpcSize);

         SPEEX_COPY(st-&gt;sw, st-&gt;winBuf, diff);
         SPEEX_COPY(st-&gt;sw+diff, in, st-&gt;frameSize-diff);
         filter_mem16(st-&gt;sw, bw_lpc1, bw_lpc2, st-&gt;sw, st-&gt;frameSize, st-&gt;lpcSize, st-&gt;mem_sw_whole, stack);

         open_loop_nbest_pitch(st-&gt;sw, st-&gt;min_pitch, st-&gt;max_pitch, st-&gt;frameSize, 
                               nol_pitch, nol_pitch_coef, 6, stack);
         ol_pitch=nol_pitch[0];
         ol_pitch_coef = nol_pitch_coef[0];
         /*Try to remove pitch multiples*/
         for (i=1;i&lt;6;i++)
         {
#ifdef FIXED_POINT
            if ((nol_pitch_coef[i]&gt;MULT16_16_Q15(nol_pitch_coef[0],27853)) &amp;&amp; 
#else
            if ((nol_pitch_coef[i]&gt;.85*nol_pitch_coef[0]) &amp;&amp; 
#endif
                (ABS(2*nol_pitch[i]-ol_pitch)&lt;=2 || ABS(3*nol_pitch[i]-ol_pitch)&lt;=3 || 
                 ABS(4*nol_pitch[i]-ol_pitch)&lt;=4 || ABS(5*nol_pitch[i]-ol_pitch)&lt;=5))
            {
               /*ol_pitch_coef=nol_pitch_coef[i];*/
               ol_pitch = nol_pitch[i];
            }
         }
         /*if (ol_pitch&gt;50)
           ol_pitch/=2;*/
         /*ol_pitch_coef = sqrt(ol_pitch_coef);*/

      } else {
         ol_pitch=0;
         ol_pitch_coef=0;
      }
      
      /*Compute &quot;real&quot; excitation*/
      SPEEX_COPY(st-&gt;exc, st-&gt;winBuf, diff);
      SPEEX_COPY(st-&gt;exc+diff, in, st-&gt;frameSize-diff);
      fir_mem16(st-&gt;exc, interp_lpc, st-&gt;exc, st-&gt;frameSize, st-&gt;lpcSize, st-&gt;mem_exc, stack);

      /* Compute open-loop excitation gain */
      {
         spx_word16_t g = compute_rms16(st-&gt;exc, st-&gt;frameSize);
         if (st-&gt;submodeID!=1 &amp;&amp; ol_pitch&gt;0)
            ol_gain = MULT16_16(g, MULT16_16_Q14(QCONST16(1.1,14),
                                spx_sqrt(QCONST32(1.,28)-MULT16_32_Q15(QCONST16(.8,15),SHL32(MULT16_16(ol_pitch_coef,ol_pitch_coef),16)))));
         else
            ol_gain = SHL32(EXTEND32(g),SIG_SHIFT);
      }
   }

#ifdef VORBIS_PSYCHO
   SPEEX_MOVE(st-&gt;psy_window, st-&gt;psy_window+st-&gt;frameSize, 256-st-&gt;frameSize);
   SPEEX_COPY(&amp;st-&gt;psy_window[256-st-&gt;frameSize], in, st-&gt;frameSize);
   compute_curve(st-&gt;psy, st-&gt;psy_window, st-&gt;curve);
   /*print_vec(st-&gt;curve, 128, &quot;curve&quot;);*/
   if (st-&gt;first)
      SPEEX_COPY(st-&gt;old_curve, st-&gt;curve, 128);
#endif

   /*VBR stuff*/
#ifndef DISABLE_VBR
   if (st-&gt;vbr &amp;&amp; (st-&gt;vbr_enabled||st-&gt;vad_enabled))
   {
      float lsp_dist=0;
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         lsp_dist += (st-&gt;old_lsp[i] - lsp[i])*(st-&gt;old_lsp[i] - lsp[i]);
      lsp_dist /= LSP_SCALING*LSP_SCALING;
      
      if (st-&gt;abr_enabled)
      {
         float qual_change=0;
         if (st-&gt;abr_drift2 * st-&gt;abr_drift &gt; 0)
         {
            /* Only adapt if long-term and short-term drift are the same sign */
            qual_change = -.00001*st-&gt;abr_drift/(1+st-&gt;abr_count);
            if (qual_change&gt;.05)
               qual_change=.05;
            if (qual_change&lt;-.05)
               qual_change=-.05;
         }
         st-&gt;vbr_quality += qual_change;
         if (st-&gt;vbr_quality&gt;10)
            st-&gt;vbr_quality=10;
         if (st-&gt;vbr_quality&lt;0)
            st-&gt;vbr_quality=0;
      }

      st-&gt;relative_quality = vbr_analysis(st-&gt;vbr, in, st-&gt;frameSize, ol_pitch, GAIN_SCALING_1*ol_pitch_coef);
      /*if (delta_qual&lt;0)*/
      /*  delta_qual*=.1*(3+st-&gt;vbr_quality);*/
      if (st-&gt;vbr_enabled) 
      {
         spx_int32_t mode;
         int choice=0;
         float min_diff=100;
         mode = 8;
         while (mode)
         {
            int v1;
            float thresh;
            v1=(int)floor(st-&gt;vbr_quality);
            if (v1==10)
               thresh = vbr_nb_thresh[mode][v1];
            else
               thresh = (st-&gt;vbr_quality-v1)*vbr_nb_thresh[mode][v1+1] + (1+v1-st-&gt;vbr_quality)*vbr_nb_thresh[mode][v1];
            if (st-&gt;relative_quality &gt; thresh &amp;&amp; 
                st-&gt;relative_quality-thresh&lt;min_diff)
            {
               choice = mode;
               min_diff = st-&gt;relative_quality-thresh;
            }
            mode--;
         }
         mode=choice;
         if (mode==0)
         {
            if (st-&gt;dtx_count==0 || lsp_dist&gt;.05 || !st-&gt;dtx_enabled || st-&gt;dtx_count&gt;20)
            {
               mode=1;
               st-&gt;dtx_count=1;
            } else {
               mode=0;
               st-&gt;dtx_count++;
            }
         } else {
            st-&gt;dtx_count=0;
         }

         speex_encoder_ctl(state, SPEEX_SET_MODE, &amp;mode);
         if (st-&gt;vbr_max&gt;0)
         {
            spx_int32_t rate;
            speex_encoder_ctl(state, SPEEX_GET_BITRATE, &amp;rate);
            if (rate &gt; st-&gt;vbr_max)
            {
               rate = st-&gt;vbr_max;
               speex_encoder_ctl(state, SPEEX_SET_BITRATE, &amp;rate);
            }
         }
         
         if (st-&gt;abr_enabled)
         {
            spx_int32_t bitrate;
            speex_encoder_ctl(state, SPEEX_GET_BITRATE, &amp;bitrate);
            st-&gt;abr_drift+=(bitrate-st-&gt;abr_enabled);
            st-&gt;abr_drift2 = .95*st-&gt;abr_drift2 + .05*(bitrate-st-&gt;abr_enabled);
            st-&gt;abr_count += 1.0;
         }

      } else {
         /*VAD only case*/
         int mode;
         if (st-&gt;relative_quality&lt;2)
         {
            if (st-&gt;dtx_count==0 || lsp_dist&gt;.05 || !st-&gt;dtx_enabled || st-&gt;dtx_count&gt;20)
            {
               st-&gt;dtx_count=1;
               mode=1;
            } else {
               mode=0;
               st-&gt;dtx_count++;
            }
         } else {
            st-&gt;dtx_count = 0;
            mode=st-&gt;submodeSelect;
         }
         /*speex_encoder_ctl(state, SPEEX_SET_MODE, &amp;mode);*/
         st-&gt;submodeID=mode;
      } 
   } else {
      st-&gt;relative_quality = -1;
   }
#endif /* #ifndef DISABLE_VBR */

   if (st-&gt;encode_submode)
   {
      /* First, transmit a zero for narrowband */
      speex_bits_pack(bits, 0, 1);

      /* Transmit the sub-mode we use for this frame */
      speex_bits_pack(bits, st-&gt;submodeID, NB_SUBMODE_BITS);

   }

   /* If null mode (no transmission), just set a couple things to zero*/
   if (st-&gt;submodes[st-&gt;submodeID] == NULL)
   {
      for (i=0;i&lt;st-&gt;frameSize;i++)
         st-&gt;exc[i]=st-&gt;sw[i]=VERY_SMALL;

      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;mem_sw[i]=0;
      st-&gt;first=1;
      st-&gt;bounded_pitch = 1;

      SPEEX_COPY(st-&gt;winBuf, in+2*st-&gt;frameSize-st-&gt;windowSize, st-&gt;windowSize-st-&gt;frameSize);

      /* Clear memory (no need to really compute it) */
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;mem_sp[i] = 0;
      return 0;

   }

   /* LSP Quantization */
   if (st-&gt;first)
   {
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;old_lsp[i] = lsp[i];
   }


   /*Quantize LSPs*/
#if 1 /*0 for unquantized*/
   SUBMODE(lsp_quant)(lsp, qlsp, st-&gt;lpcSize, bits);
#else
   for (i=0;i&lt;st-&gt;lpcSize;i++)
     qlsp[i]=lsp[i];
#endif

   /*If we use low bit-rate pitch mode, transmit open-loop pitch*/
   if (SUBMODE(lbr_pitch)!=-1)
   {
      speex_bits_pack(bits, ol_pitch-st-&gt;min_pitch, 7);
   } 

   if (SUBMODE(forced_pitch_gain))
   {
      int quant;
      /* This just damps the pitch a bit, because it tends to be too aggressive when forced */
      ol_pitch_coef = MULT16_16_Q15(QCONST16(.9,15), ol_pitch_coef);
#ifdef FIXED_POINT
      quant = PSHR16(MULT16_16_16(15, ol_pitch_coef),GAIN_SHIFT);
#else
      quant = (int)floor(.5+15*ol_pitch_coef*GAIN_SCALING_1);
#endif
      if (quant&gt;15)
         quant=15;
      if (quant&lt;0)
         quant=0;
      speex_bits_pack(bits, quant, 4);
      ol_pitch_coef=MULT16_16_P15(QCONST16(0.066667,15),SHL16(quant,GAIN_SHIFT));
   }
   
   
   /*Quantize and transmit open-loop excitation gain*/
#ifdef FIXED_POINT
   {
      int qe = scal_quant32(ol_gain, ol_gain_table, 32);
      /*ol_gain = exp(qe/3.5)*SIG_SCALING;*/
      ol_gain = MULT16_32_Q15(28406,ol_gain_table[qe]);
      speex_bits_pack(bits, qe, 5);
   }
#else
   {
      int qe = (int)(floor(.5+3.5*log(ol_gain*1.0/SIG_SCALING)));
      if (qe&lt;0)
         qe=0;
      if (qe&gt;31)
         qe=31;
      ol_gain = exp(qe/3.5)*SIG_SCALING;
      speex_bits_pack(bits, qe, 5);
   }
#endif



   /* Special case for first frame */
   if (st-&gt;first)
   {
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;old_qlsp[i] = qlsp[i];
   }

   /* Target signal */
   ALLOC(target, st-&gt;subframeSize, spx_word16_t);
   ALLOC(innov, st-&gt;subframeSize, spx_sig_t);
   ALLOC(exc32, st-&gt;subframeSize, spx_word32_t);
   ALLOC(ringing, st-&gt;subframeSize, spx_word16_t);
   ALLOC(syn_resp, st-&gt;subframeSize, spx_word16_t);
   ALLOC(real_exc, st-&gt;subframeSize, spx_word16_t);
   ALLOC(mem, st-&gt;lpcSize, spx_mem_t);

   /* Loop on sub-frames */
   for (sub=0;sub&lt;st-&gt;nbSubframes;sub++)
   {
      int   offset;
      spx_word16_t *sw;
      spx_word16_t *exc;
      int pitch;
      int response_bound = st-&gt;subframeSize;

      /* Offset relative to start of frame */
      offset = st-&gt;subframeSize*sub;
      /* Excitation */
      exc=st-&gt;exc+offset;
      /* Weighted signal */
      sw=st-&gt;sw+offset;
      
      /* LSP interpolation (quantized and unquantized) */
      lsp_interpolate(st-&gt;old_lsp, lsp, interp_lsp, st-&gt;lpcSize, sub, st-&gt;nbSubframes);
      lsp_interpolate(st-&gt;old_qlsp, qlsp, interp_qlsp, st-&gt;lpcSize, sub, st-&gt;nbSubframes);

      /* Make sure the filters are stable */
      lsp_enforce_margin(interp_lsp, st-&gt;lpcSize, LSP_MARGIN);
      lsp_enforce_margin(interp_qlsp, st-&gt;lpcSize, LSP_MARGIN);

      /* Compute interpolated LPCs (quantized and unquantized) */
      lsp_to_lpc(interp_lsp, interp_lpc, st-&gt;lpcSize,stack);

      lsp_to_lpc(interp_qlsp, interp_qlpc, st-&gt;lpcSize, stack);

      /* Compute analysis filter gain at w=pi (for use in SB-CELP) */
      {
         spx_word32_t pi_g=LPC_SCALING;
         for (i=0;i&lt;st-&gt;lpcSize;i+=2)
         {
            /*pi_g += -st-&gt;interp_qlpc[i] +  st-&gt;interp_qlpc[i+1];*/
            pi_g = ADD32(pi_g, SUB32(EXTEND32(interp_qlpc[i+1]),EXTEND32(interp_qlpc[i])));
         }
         st-&gt;pi_gain[sub] = pi_g;
      }

#ifdef VORBIS_PSYCHO
      {
         float curr_curve[128];
         float fact = ((float)sub+1.0f)/st-&gt;nbSubframes;
         for (i=0;i&lt;128;i++)
            curr_curve[i] = (1.0f-fact)*st-&gt;old_curve[i] + fact*st-&gt;curve[i];
         curve_to_lpc(st-&gt;psy, curr_curve, bw_lpc1, bw_lpc2, 10);
      }
#else
      /* Compute bandwidth-expanded (unquantized) LPCs for perceptual weighting */
      bw_lpc(st-&gt;gamma1, interp_lpc, bw_lpc1, st-&gt;lpcSize);
      if (st-&gt;gamma2&gt;=0)
         bw_lpc(st-&gt;gamma2, interp_lpc, bw_lpc2, st-&gt;lpcSize);
      else
      {
         for (i=0;i&lt;st-&gt;lpcSize;i++)
            bw_lpc2[i]=0;
      }
      /*print_vec(st-&gt;bw_lpc1, 10, &quot;bw_lpc&quot;);*/
#endif

      /*FIXME: This will break if we change the window size */
      speex_assert(st-&gt;windowSize-st-&gt;frameSize == st-&gt;subframeSize);
      if (sub==0)
      {
         for (i=0;i&lt;st-&gt;subframeSize;i++)
            real_exc[i] = sw[i] = st-&gt;winBuf[i];
      } else {
         for (i=0;i&lt;st-&gt;subframeSize;i++)
            real_exc[i] = sw[i] = in[i+((sub-1)*st-&gt;subframeSize)];
      }
      fir_mem16(real_exc, interp_qlpc, real_exc, st-&gt;subframeSize, st-&gt;lpcSize, st-&gt;mem_exc2, stack);
      
      if (st-&gt;complexity==0)
         response_bound &gt;&gt;= 1;
      compute_impulse_response(interp_qlpc, bw_lpc1, bw_lpc2, syn_resp, response_bound, st-&gt;lpcSize, stack);
      for (i=response_bound;i&lt;st-&gt;subframeSize;i++)
         syn_resp[i]=VERY_SMALL;
      
      /* Compute zero response of A(z/g1) / ( A(z/g2) * A(z) ) */
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         mem[i]=SHL32(st-&gt;mem_sp[i],1);
      for (i=0;i&lt;st-&gt;subframeSize;i++)
         ringing[i] = VERY_SMALL;
#ifdef SHORTCUTS2
      iir_mem16(ringing, interp_qlpc, ringing, response_bound, st-&gt;lpcSize, mem, stack);
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         mem[i]=SHL32(st-&gt;mem_sw[i],1);
      filter_mem16(ringing, st-&gt;bw_lpc1, st-&gt;bw_lpc2, ringing, response_bound, st-&gt;lpcSize, mem, stack);
      SPEEX_MEMSET(&amp;ringing[response_bound], 0, st-&gt;subframeSize-response_bound);
#else
      iir_mem16(ringing, interp_qlpc, ringing, st-&gt;subframeSize, st-&gt;lpcSize, mem, stack);
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         mem[i]=SHL32(st-&gt;mem_sw[i],1);
      filter_mem16(ringing, bw_lpc1, bw_lpc2, ringing, st-&gt;subframeSize, st-&gt;lpcSize, mem, stack);
#endif
      
      /* Compute weighted signal */
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         mem[i]=st-&gt;mem_sw[i];
      filter_mem16(sw, bw_lpc1, bw_lpc2, sw, st-&gt;subframeSize, st-&gt;lpcSize, mem, stack);
      
      if (st-&gt;complexity==0)
         for (i=0;i&lt;st-&gt;lpcSize;i++)
            st-&gt;mem_sw[i]=mem[i];
      
      /* Compute target signal (saturation prevents overflows on clipped input speech) */
      for (i=0;i&lt;st-&gt;subframeSize;i++)
         target[i]=EXTRACT16(SATURATE(SUB32(sw[i],PSHR32(ringing[i],1)),32767));

      /* Reset excitation */
      SPEEX_MEMSET(exc, 0, st-&gt;subframeSize);

      /* If we have a long-term predictor (otherwise, something&#39;s wrong) */
      speex_assert (SUBMODE(ltp_quant));
      {
         int pit_min, pit_max;
         /* Long-term prediction */
         if (SUBMODE(lbr_pitch) != -1)
         {
            /* Low bit-rate pitch handling */
            int margin;
            margin = SUBMODE(lbr_pitch);
            if (margin)
            {
               if (ol_pitch &lt; st-&gt;min_pitch+margin-1)
                  ol_pitch=st-&gt;min_pitch+margin-1;
               if (ol_pitch &gt; st-&gt;max_pitch-margin)
                  ol_pitch=st-&gt;max_pitch-margin;
               pit_min = ol_pitch-margin+1;
               pit_max = ol_pitch+margin;
            } else {
               pit_min=pit_max=ol_pitch;
            }
         } else {
            pit_min = st-&gt;min_pitch;
            pit_max = st-&gt;max_pitch;
         }
         
         /* Force pitch to use only the current frame if needed */
         if (st-&gt;bounded_pitch &amp;&amp; pit_max&gt;offset)
            pit_max=offset;

         /* Perform pitch search */
         pitch = SUBMODE(ltp_quant)(target, sw, interp_qlpc, bw_lpc1, bw_lpc2,
                                    exc32, SUBMODE(ltp_params), pit_min, pit_max, ol_pitch_coef,
                                    st-&gt;lpcSize, st-&gt;subframeSize, bits, stack, 
                                    exc, syn_resp, st-&gt;complexity, 0, st-&gt;plc_tuning, &amp;st-&gt;cumul_gain);

         st-&gt;pitch[sub]=pitch;
      }
      /* Quantization of innovation */
      SPEEX_MEMSET(innov, 0, st-&gt;subframeSize);
      
      /* FIXME: Make sure this is save from overflows (so far so good) */
      for (i=0;i&lt;st-&gt;subframeSize;i++)
         real_exc[i] = EXTRACT16(SUB32(EXTEND32(real_exc[i]), PSHR32(exc32[i],SIG_SHIFT-1)));
      
      ener = SHL32(EXTEND32(compute_rms16(real_exc, st-&gt;subframeSize)),SIG_SHIFT);
      
      /*FIXME: Should use DIV32_16 and make sure result fits in 16 bits */
#ifdef FIXED_POINT
      {
         spx_word32_t f = PDIV32(ener,PSHR32(ol_gain,SIG_SHIFT));
         if (f&lt;=32767)
            fine_gain = f;
         else
            fine_gain = 32767;
      }
#else
      fine_gain = PDIV32_16(ener,PSHR32(ol_gain,SIG_SHIFT));
#endif
      /* Calculate gain correction for the sub-frame (if any) */
      if (SUBMODE(have_subframe_gain)) 
      {
         int qe;
         if (SUBMODE(have_subframe_gain)==3)
         {
            qe = scal_quant(fine_gain, exc_gain_quant_scal3_bound, 8);
            speex_bits_pack(bits, qe, 3);
            ener=MULT16_32_Q14(exc_gain_quant_scal3[qe],ol_gain);
         } else {
            qe = scal_quant(fine_gain, exc_gain_quant_scal1_bound, 2);
            speex_bits_pack(bits, qe, 1);
            ener=MULT16_32_Q14(exc_gain_quant_scal1[qe],ol_gain);               
         }
      } else {
         ener=ol_gain;
      }
      
      /*printf (&quot;%f %f\n&quot;, ener, ol_gain);*/
      
      /* Normalize innovation */
      signal_div(target, target, ener, st-&gt;subframeSize);
      
      /* Quantize innovation */
      speex_assert (SUBMODE(innovation_quant));
      {
         /* Codebook search */
         SUBMODE(innovation_quant)(target, interp_qlpc, bw_lpc1, bw_lpc2, 
                  SUBMODE(innovation_params), st-&gt;lpcSize, st-&gt;subframeSize, 
                  innov, syn_resp, bits, stack, st-&gt;complexity, SUBMODE(double_codebook));
         
         /* De-normalize innovation and update excitation */
         signal_mul(innov, innov, ener, st-&gt;subframeSize);
         
         for (i=0;i&lt;st-&gt;subframeSize;i++)
            exc[i] = EXTRACT16(SATURATE32(PSHR32(ADD32(SHL32(exc32[i],1),innov[i]),SIG_SHIFT),32767));

         /* In some (rare) modes, we do a second search (more bits) to reduce noise even more */
         if (SUBMODE(double_codebook)) {
            char *tmp_stack=stack;
            VARDECL(spx_sig_t *innov2);
            ALLOC(innov2, st-&gt;subframeSize, spx_sig_t);
            SPEEX_MEMSET(innov2, 0, st-&gt;subframeSize);
            for (i=0;i&lt;st-&gt;subframeSize;i++)
               target[i]=MULT16_16_P13(QCONST16(2.2f,13), target[i]);
            SUBMODE(innovation_quant)(target, interp_qlpc, bw_lpc1, bw_lpc2, 
                                      SUBMODE(innovation_params), st-&gt;lpcSize, st-&gt;subframeSize, 
                                      innov2, syn_resp, bits, stack, st-&gt;complexity, 0);
            signal_mul(innov2, innov2, MULT16_32_Q15(QCONST16(0.454545f,15),ener), st-&gt;subframeSize);
            for (i=0;i&lt;st-&gt;subframeSize;i++)
               innov[i] = ADD32(innov[i],innov2[i]);
            stack = tmp_stack;
         }
         for (i=0;i&lt;st-&gt;subframeSize;i++)
            exc[i] = EXTRACT16(SATURATE32(PSHR32(ADD32(SHL32(exc32[i],1),innov[i]),SIG_SHIFT),32767));
         if (st-&gt;innov_rms_save)
         {
            st-&gt;innov_rms_save[sub] = compute_rms(innov, st-&gt;subframeSize);
         }
      }

      /* Final signal synthesis from excitation */
      iir_mem16(exc, interp_qlpc, sw, st-&gt;subframeSize, st-&gt;lpcSize, st-&gt;mem_sp, stack);

      /* Compute weighted signal again, from synthesized speech (not sure it&#39;s the right thing) */
      if (st-&gt;complexity!=0)
         filter_mem16(sw, bw_lpc1, bw_lpc2, sw, st-&gt;subframeSize, st-&gt;lpcSize, st-&gt;mem_sw, stack);
      
   }

   /* Store the LSPs for interpolation in the next frame */
   if (st-&gt;submodeID&gt;=1)
   {
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;old_lsp[i] = lsp[i];
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;old_qlsp[i] = qlsp[i];
   }

#ifdef VORBIS_PSYCHO
   if (st-&gt;submodeID&gt;=1)
      SPEEX_COPY(st-&gt;old_curve, st-&gt;curve, 128);
#endif

   if (st-&gt;submodeID==1)
   {
#ifndef DISABLE_VBR
      if (st-&gt;dtx_count)
         speex_bits_pack(bits, 15, 4);
      else
#endif
         speex_bits_pack(bits, 0, 4);
   }

   /* The next frame will not be the first (Duh!) */
   st-&gt;first = 0;
   SPEEX_COPY(st-&gt;winBuf, in+2*st-&gt;frameSize-st-&gt;windowSize, st-&gt;windowSize-st-&gt;frameSize);

   if (SUBMODE(innovation_quant) == noise_codebook_quant || st-&gt;submodeID==0)
      st-&gt;bounded_pitch = 1;
   else
      st-&gt;bounded_pitch = 0;

   return 1;
}

void *nb_decoder_init(const SpeexMode *m)
{
   DecState *st;
   const SpeexNBMode *mode;
   int i;

   mode=(const SpeexNBMode*)m-&gt;mode;
   st = (DecState *)speex_alloc(sizeof(DecState));
   if (!st)
      return NULL;
#if defined(VAR_ARRAYS) || defined (USE_ALLOCA)
   st-&gt;stack = NULL;
#else
   st-&gt;stack = (char*)speex_alloc_scratch(NB_DEC_STACK);
#endif

   st-&gt;mode=m;


   st-&gt;encode_submode = 1;

   st-&gt;first=1;
   /* Codec parameters, should eventually have several &quot;modes&quot;*/
   st-&gt;frameSize = mode-&gt;frameSize;
   st-&gt;nbSubframes=mode-&gt;frameSize/mode-&gt;subframeSize;
   st-&gt;subframeSize=mode-&gt;subframeSize;
   st-&gt;lpcSize = mode-&gt;lpcSize;
   st-&gt;min_pitch=mode-&gt;pitchStart;
   st-&gt;max_pitch=mode-&gt;pitchEnd;

   st-&gt;submodes=mode-&gt;submodes;
   st-&gt;submodeID=mode-&gt;defaultSubmode;

   st-&gt;lpc_enh_enabled=1;

   st-&gt;excBuf = (spx_word16_t*)speex_alloc((st-&gt;frameSize + 2*st-&gt;max_pitch + st-&gt;subframeSize + 12)*sizeof(spx_word16_t));
   st-&gt;exc = st-&gt;excBuf + 2*st-&gt;max_pitch + st-&gt;subframeSize + 6;
   SPEEX_MEMSET(st-&gt;excBuf, 0, st-&gt;frameSize + st-&gt;max_pitch);

   st-&gt;interp_qlpc = (spx_coef_t*)speex_alloc(st-&gt;lpcSize*sizeof(spx_coef_t));
   st-&gt;old_qlsp = (spx_lsp_t*)speex_alloc(st-&gt;lpcSize*sizeof(spx_lsp_t));
   st-&gt;mem_sp = (spx_mem_t*)speex_alloc(st-&gt;lpcSize*sizeof(spx_mem_t));
   st-&gt;pi_gain = (spx_word32_t*)speex_alloc((st-&gt;nbSubframes)*sizeof(spx_word32_t));
   st-&gt;last_pitch = 40;
   st-&gt;count_lost=0;
   st-&gt;pitch_gain_buf[0] = st-&gt;pitch_gain_buf[1] = st-&gt;pitch_gain_buf[2] = 0;
   st-&gt;pitch_gain_buf_idx = 0;
   st-&gt;seed = 1000;
   
   st-&gt;sampling_rate=8000;
   st-&gt;last_ol_gain = 0;

   st-&gt;user_callback.func = &amp;speex_default_user_handler;
   st-&gt;user_callback.data = NULL;
   for (i=0;i&lt;16;i++)
      st-&gt;speex_callbacks[i].func = NULL;

   st-&gt;voc_m1=st-&gt;voc_m2=st-&gt;voc_mean=0;
   st-&gt;voc_offset=0;
   st-&gt;dtx_enabled=0;
   st-&gt;isWideband = 0;
   st-&gt;highpass_enabled = 1;

#ifdef ENABLE_VALGRIND
   VALGRIND_MAKE_READABLE(st, NB_DEC_STACK);
#endif
   return st;
}

void nb_decoder_destroy(void *state)
{
   DecState *st;
   st=(DecState*)state;
   
#if !(defined(VAR_ARRAYS) || defined (USE_ALLOCA))
   speex_free_scratch(st-&gt;stack);
#endif

   speex_free (st-&gt;excBuf);
   speex_free (st-&gt;interp_qlpc);
   speex_free (st-&gt;old_qlsp);
   speex_free (st-&gt;mem_sp);
   speex_free (st-&gt;pi_gain);

   speex_free(state);
}

#define median3(a, b, c)	((a) &lt; (b) ? ((b) &lt; (c) ? (b) : ((a) &lt; (c) ? (c) : (a))) : ((c) &lt; (b) ? (b) : ((c) &lt; (a) ? (c) : (a))))

#ifdef FIXED_POINT
const spx_word16_t attenuation[10] = {32767, 31483, 27923, 22861, 17278, 12055, 7764, 4616, 2533, 1283};
#else
const spx_word16_t attenuation[10] = {1., 0.961, 0.852, 0.698, 0.527, 0.368, 0.237, 0.141, 0.077, 0.039};

#endif

static void nb_decode_lost(DecState *st, spx_word16_t *out, char *stack)
{
   int i;
   int pitch_val;
   spx_word16_t pitch_gain;
   spx_word16_t fact;
   spx_word16_t gain_med;
   spx_word16_t innov_gain;
   spx_word16_t noise_gain;
   
   if (st-&gt;count_lost&lt;10)
      fact = attenuation[st-&gt;count_lost];
   else
      fact = 0;

   gain_med = median3(st-&gt;pitch_gain_buf[0], st-&gt;pitch_gain_buf[1], st-&gt;pitch_gain_buf[2]);
   if (gain_med &lt; st-&gt;last_pitch_gain)
      st-&gt;last_pitch_gain = gain_med;
   
#ifdef FIXED_POINT
   pitch_gain = st-&gt;last_pitch_gain;
   if (pitch_gain&gt;54)
      pitch_gain = 54;
   pitch_gain = SHL16(pitch_gain, 9);
#else   
   pitch_gain = GAIN_SCALING_1*st-&gt;last_pitch_gain;
   if (pitch_gain&gt;.85)
      pitch_gain=.85;
#endif
   pitch_gain = MULT16_16_Q15(fact,pitch_gain) + VERY_SMALL;
   /* FIXME: This was rms of innovation (not exc) */
   innov_gain = compute_rms16(st-&gt;exc, st-&gt;frameSize);
   noise_gain = MULT16_16_Q15(innov_gain, MULT16_16_Q15(fact, SUB16(Q15ONE,MULT16_16_Q15(pitch_gain,pitch_gain))));
   /* Shift all buffers by one frame */
   SPEEX_MOVE(st-&gt;excBuf, st-&gt;excBuf+st-&gt;frameSize, 2*st-&gt;max_pitch + st-&gt;subframeSize + 12);
   

   pitch_val = st-&gt;last_pitch + SHR32((spx_int32_t)speex_rand(1+st-&gt;count_lost, &amp;st-&gt;seed),SIG_SHIFT);
   if (pitch_val &gt; st-&gt;max_pitch)
      pitch_val = st-&gt;max_pitch;
   if (pitch_val &lt; st-&gt;min_pitch)
      pitch_val = st-&gt;min_pitch;
   for (i=0;i&lt;st-&gt;frameSize;i++)
   {
      st-&gt;exc[i]= MULT16_16_Q15(pitch_gain, (st-&gt;exc[i-pitch_val]+VERY_SMALL)) + 
            speex_rand(noise_gain, &amp;st-&gt;seed);
   }

   bw_lpc(QCONST16(.98,15), st-&gt;interp_qlpc, st-&gt;interp_qlpc, st-&gt;lpcSize);
   iir_mem16(&amp;st-&gt;exc[-st-&gt;subframeSize], st-&gt;interp_qlpc, out, st-&gt;frameSize,
             st-&gt;lpcSize, st-&gt;mem_sp, stack);
   highpass(out, out, st-&gt;frameSize, HIGHPASS_NARROWBAND|HIGHPASS_OUTPUT, st-&gt;mem_hp);
   
   st-&gt;first = 0;
   st-&gt;count_lost++;
   st-&gt;pitch_gain_buf[st-&gt;pitch_gain_buf_idx++] = PSHR16(pitch_gain,9);
   if (st-&gt;pitch_gain_buf_idx &gt; 2) /* rollover */
      st-&gt;pitch_gain_buf_idx = 0;
}

/* Just so we don&#39;t need to carry the complete wideband mode information */
static const int wb_skip_table[8] = {0, 36, 112, 192, 352, 0, 0, 0};
   
int nb_decode(void *state, SpeexBits *bits, void *vout)
{
   DecState *st;
   int i, sub;
   int pitch;
   spx_word16_t pitch_gain[3];
   spx_word32_t ol_gain=0;
   int ol_pitch=0;
   spx_word16_t ol_pitch_coef=0;
   int best_pitch=40;
   spx_word16_t best_pitch_gain=0;
   int wideband;
   int m;
   char *stack;
   VARDECL(spx_sig_t *innov);
   VARDECL(spx_word32_t *exc32);
   VARDECL(spx_coef_t *ak);
   VARDECL(spx_lsp_t *qlsp);
   spx_word16_t pitch_average=0;
   
   spx_word16_t *out = (spx_word16_t*)vout;
   VARDECL(spx_lsp_t *interp_qlsp);

   st=(DecState*)state;
   stack=st-&gt;stack;

   /* Check if we&#39;re in DTX mode*/
   if (!bits &amp;&amp; st-&gt;dtx_enabled)
   {
      st-&gt;submodeID=0;
   } else 
   {
      /* If bits is NULL, consider the packet to be lost (what could we do anyway) */
      if (!bits)
      {
         nb_decode_lost(st, out, stack);
         return 0;
      }

      if (st-&gt;encode_submode)
      {

      /* Search for next narrowband block (handle requests, skip wideband blocks) */
      do {
         if (speex_bits_remaining(bits)&lt;5)
            return -1;
         wideband = speex_bits_unpack_unsigned(bits, 1);
         if (wideband) /* Skip wideband block (for compatibility) */
         {
            int submode;
            int advance;
            advance = submode = speex_bits_unpack_unsigned(bits, SB_SUBMODE_BITS);
            /*speex_mode_query(&amp;speex_wb_mode, SPEEX_SUBMODE_BITS_PER_FRAME, &amp;advance);*/
            advance = wb_skip_table[submode];
            if (advance &lt; 0)
            {
               speex_notify(&quot;Invalid mode encountered. The stream is corrupted.&quot;);
               return -2;
            } 
            advance -= (SB_SUBMODE_BITS+1);
            speex_bits_advance(bits, advance);
            
            if (speex_bits_remaining(bits)&lt;5)
               return -1;
            wideband = speex_bits_unpack_unsigned(bits, 1);
            if (wideband)
            {
               advance = submode = speex_bits_unpack_unsigned(bits, SB_SUBMODE_BITS);
               /*speex_mode_query(&amp;speex_wb_mode, SPEEX_SUBMODE_BITS_PER_FRAME, &amp;advance);*/
               advance = wb_skip_table[submode];
               if (advance &lt; 0)
               {
                  speex_notify(&quot;Invalid mode encountered. The stream is corrupted.&quot;);
                  return -2;
               } 
               advance -= (SB_SUBMODE_BITS+1);
               speex_bits_advance(bits, advance);
               wideband = speex_bits_unpack_unsigned(bits, 1);
               if (wideband)
               {
                  speex_notify(&quot;More than two wideband layers found. The stream is corrupted.&quot;);
                  return -2;
               }

            }
         }
         if (speex_bits_remaining(bits)&lt;4)
            return -1;
         /* FIXME: Check for overflow */
         m = speex_bits_unpack_unsigned(bits, 4);
         if (m==15) /* We found a terminator */
         {
            return -1;
         } else if (m==14) /* Speex in-band request */
         {
            int ret = speex_inband_handler(bits, st-&gt;speex_callbacks, state);
            if (ret)
               return ret;
         } else if (m==13) /* User in-band request */
         {
            int ret = st-&gt;user_callback.func(bits, state, st-&gt;user_callback.data);
            if (ret)
               return ret;
         } else if (m&gt;8) /* Invalid mode */
         {
            speex_notify(&quot;Invalid mode encountered. The stream is corrupted.&quot;);
            return -2;
         }
      
      } while (m&gt;8);

      /* Get the sub-mode that was used */
      st-&gt;submodeID = m;
      }

   }

   /* Shift all buffers by one frame */
   SPEEX_MOVE(st-&gt;excBuf, st-&gt;excBuf+st-&gt;frameSize, 2*st-&gt;max_pitch + st-&gt;subframeSize + 12);

   /* If null mode (no transmission), just set a couple things to zero*/
   if (st-&gt;submodes[st-&gt;submodeID] == NULL)
   {
      VARDECL(spx_coef_t *lpc);
      ALLOC(lpc, st-&gt;lpcSize, spx_coef_t);
      bw_lpc(QCONST16(0.93f,15), st-&gt;interp_qlpc, lpc, st-&gt;lpcSize);
      {
         spx_word16_t innov_gain=0;
         /* FIXME: This was innov, not exc */
         innov_gain = compute_rms16(st-&gt;exc, st-&gt;frameSize);
         for (i=0;i&lt;st-&gt;frameSize;i++)
            st-&gt;exc[i]=speex_rand(innov_gain, &amp;st-&gt;seed);
      }


      st-&gt;first=1;

      /* Final signal synthesis from excitation */
      iir_mem16(st-&gt;exc, lpc, out, st-&gt;frameSize, st-&gt;lpcSize, st-&gt;mem_sp, stack);

      st-&gt;count_lost=0;
      return 0;
   }

   ALLOC(qlsp, st-&gt;lpcSize, spx_lsp_t);

   /* Unquantize LSPs */
   SUBMODE(lsp_unquant)(qlsp, st-&gt;lpcSize, bits);

   /*Damp memory if a frame was lost and the LSP changed too much*/
   if (st-&gt;count_lost)
   {
      spx_word16_t fact;
      spx_word32_t lsp_dist=0;
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         lsp_dist = ADD32(lsp_dist, EXTEND32(ABS(st-&gt;old_qlsp[i] - qlsp[i])));
#ifdef FIXED_POINT
      fact = SHR16(19661,SHR32(lsp_dist,LSP_SHIFT+2));      
#else
      fact = .6*exp(-.2*lsp_dist);
#endif
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;mem_sp[i] = MULT16_32_Q15(fact,st-&gt;mem_sp[i]);
   }


   /* Handle first frame and lost-packet case */
   if (st-&gt;first || st-&gt;count_lost)
   {
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;old_qlsp[i] = qlsp[i];
   }

   /* Get open-loop pitch estimation for low bit-rate pitch coding */
   if (SUBMODE(lbr_pitch)!=-1)
   {
      ol_pitch = st-&gt;min_pitch+speex_bits_unpack_unsigned(bits, 7);
   } 
   
   if (SUBMODE(forced_pitch_gain))
   {
      int quant;
      quant = speex_bits_unpack_unsigned(bits, 4);
      ol_pitch_coef=MULT16_16_P15(QCONST16(0.066667,15),SHL16(quant,GAIN_SHIFT));
   }
   
   /* Get global excitation gain */
   {
      int qe;
      qe = speex_bits_unpack_unsigned(bits, 5);
#ifdef FIXED_POINT
      /* FIXME: Perhaps we could slightly lower the gain here when the output is going to saturate? */
      ol_gain = MULT16_32_Q15(28406,ol_gain_table[qe]);
#else
      ol_gain = SIG_SCALING*exp(qe/3.5);
#endif
   }

   ALLOC(ak, st-&gt;lpcSize, spx_coef_t);
   ALLOC(innov, st-&gt;subframeSize, spx_sig_t);
   ALLOC(exc32, st-&gt;subframeSize, spx_word32_t);

   if (st-&gt;submodeID==1)
   {
      int extra;
      extra = speex_bits_unpack_unsigned(bits, 4);

      if (extra==15)
         st-&gt;dtx_enabled=1;
      else
         st-&gt;dtx_enabled=0;
   }
   if (st-&gt;submodeID&gt;1)
      st-&gt;dtx_enabled=0;

   /*Loop on subframes */
   for (sub=0;sub&lt;st-&gt;nbSubframes;sub++)
   {
      int offset;
      spx_word16_t *exc;
      spx_word16_t *sp;
      spx_word16_t *innov_save = NULL;
      spx_word16_t tmp;

      /* Offset relative to start of frame */
      offset = st-&gt;subframeSize*sub;
      /* Excitation */
      exc=st-&gt;exc+offset;
      /* Original signal */
      sp=out+offset;
      if (st-&gt;innov_save)
         innov_save = st-&gt;innov_save+offset;


      /* Reset excitation */
      SPEEX_MEMSET(exc, 0, st-&gt;subframeSize);

      /*Adaptive codebook contribution*/
      speex_assert (SUBMODE(ltp_unquant));
      {
         int pit_min, pit_max;
         /* Handle pitch constraints if any */
         if (SUBMODE(lbr_pitch) != -1)
         {
            int margin;
            margin = SUBMODE(lbr_pitch);
            if (margin)
            {
/* GT - need optimization?
               if (ol_pitch &lt; st-&gt;min_pitch+margin-1)
                  ol_pitch=st-&gt;min_pitch+margin-1;
               if (ol_pitch &gt; st-&gt;max_pitch-margin)
                  ol_pitch=st-&gt;max_pitch-margin;
               pit_min = ol_pitch-margin+1;
               pit_max = ol_pitch+margin;
*/
               pit_min = ol_pitch-margin+1;
               if (pit_min &lt; st-&gt;min_pitch)
		  pit_min = st-&gt;min_pitch;
               pit_max = ol_pitch+margin;
               if (pit_max &gt; st-&gt;max_pitch)
		  pit_max = st-&gt;max_pitch;
            } else {
               pit_min = pit_max = ol_pitch;
            }
         } else {
            pit_min = st-&gt;min_pitch;
            pit_max = st-&gt;max_pitch;
         }



         SUBMODE(ltp_unquant)(exc, exc32, pit_min, pit_max, ol_pitch_coef, SUBMODE(ltp_params), 
                 st-&gt;subframeSize, &amp;pitch, &amp;pitch_gain[0], bits, stack, 
                 st-&gt;count_lost, offset, st-&gt;last_pitch_gain, 0);

         /* Ensuring that things aren&#39;t blowing up as would happen if e.g. an encoder is 
         crafting packets to make us produce NaNs and slow down the decoder (vague DoS threat).
         We can probably be even more aggressive and limit to 15000 or so. */
         sanitize_values32(exc32, NEG32(QCONST32(32000,SIG_SHIFT-1)), QCONST32(32000,SIG_SHIFT-1), st-&gt;subframeSize);
         
         tmp = gain_3tap_to_1tap(pitch_gain);

         pitch_average += tmp;
         if ((tmp&gt;best_pitch_gain&amp;&amp;ABS(2*best_pitch-pitch)&gt;=3&amp;&amp;ABS(3*best_pitch-pitch)&gt;=4&amp;&amp;ABS(4*best_pitch-pitch)&gt;=5) 
              || (tmp&gt;MULT16_16_Q15(QCONST16(.6,15),best_pitch_gain)&amp;&amp;(ABS(best_pitch-2*pitch)&lt;3||ABS(best_pitch-3*pitch)&lt;4||ABS(best_pitch-4*pitch)&lt;5)) 
              || (MULT16_16_Q15(QCONST16(.67,15),tmp)&gt;best_pitch_gain&amp;&amp;(ABS(2*best_pitch-pitch)&lt;3||ABS(3*best_pitch-pitch)&lt;4||ABS(4*best_pitch-pitch)&lt;5)) )
         {
            best_pitch = pitch;
            if (tmp &gt; best_pitch_gain)
               best_pitch_gain = tmp;
         }
      }
      
      /* Unquantize the innovation */
      {
         int q_energy;
         spx_word32_t ener;
         
         SPEEX_MEMSET(innov, 0, st-&gt;subframeSize);

         /* Decode sub-frame gain correction */
         if (SUBMODE(have_subframe_gain)==3)
         {
            q_energy = speex_bits_unpack_unsigned(bits, 3);
            ener = MULT16_32_Q14(exc_gain_quant_scal3[q_energy],ol_gain);
         } else if (SUBMODE(have_subframe_gain)==1)
         {
            q_energy = speex_bits_unpack_unsigned(bits, 1);
            ener = MULT16_32_Q14(exc_gain_quant_scal1[q_energy],ol_gain);
         } else {
            ener = ol_gain;
         }
                  
         speex_assert (SUBMODE(innovation_unquant));
         {
            /*Fixed codebook contribution*/
            SUBMODE(innovation_unquant)(innov, SUBMODE(innovation_params), st-&gt;subframeSize, bits, stack, &amp;st-&gt;seed);
            /* De-normalize innovation and update excitation */

            signal_mul(innov, innov, ener, st-&gt;subframeSize);

            /* Decode second codebook (only for some modes) */
            if (SUBMODE(double_codebook))
            {
               char *tmp_stack=stack;
               VARDECL(spx_sig_t *innov2);
               ALLOC(innov2, st-&gt;subframeSize, spx_sig_t);
               SPEEX_MEMSET(innov2, 0, st-&gt;subframeSize);
               SUBMODE(innovation_unquant)(innov2, SUBMODE(innovation_params), st-&gt;subframeSize, bits, stack, &amp;st-&gt;seed);
               signal_mul(innov2, innov2, MULT16_32_Q15(QCONST16(0.454545f,15),ener), st-&gt;subframeSize);
               for (i=0;i&lt;st-&gt;subframeSize;i++)
                  innov[i] = ADD32(innov[i], innov2[i]);
               stack = tmp_stack;
            }
            for (i=0;i&lt;st-&gt;subframeSize;i++)
               exc[i]=EXTRACT16(SATURATE32(PSHR32(ADD32(SHL32(exc32[i],1),innov[i]),SIG_SHIFT),32767));
            /*print_vec(exc, 40, &quot;innov&quot;);*/
            if (innov_save)
            {
               for (i=0;i&lt;st-&gt;subframeSize;i++)
                  innov_save[i] = EXTRACT16(PSHR32(innov[i], SIG_SHIFT));
            }
         }

         /*Vocoder mode*/
         if (st-&gt;submodeID==1) 
         {
            spx_word16_t g=ol_pitch_coef;
            g=MULT16_16_P14(QCONST16(1.5f,14),(g-QCONST16(.2f,6)));
            if (g&lt;0)
               g=0;
            if (g&gt;GAIN_SCALING)
               g=GAIN_SCALING;
            
            SPEEX_MEMSET(exc, 0, st-&gt;subframeSize);
            while (st-&gt;voc_offset&lt;st-&gt;subframeSize)
            {
               /* exc[st-&gt;voc_offset]= g*sqrt(2*ol_pitch)*ol_gain;
                  Not quite sure why we need the factor of two in the sqrt */
               if (st-&gt;voc_offset&gt;=0)
                  exc[st-&gt;voc_offset]=MULT16_16(spx_sqrt(MULT16_16_16(2,ol_pitch)),EXTRACT16(PSHR32(MULT16_16(g,PSHR32(ol_gain,SIG_SHIFT)),6)));
               st-&gt;voc_offset+=ol_pitch;
            }
            st-&gt;voc_offset -= st-&gt;subframeSize;
            
            for (i=0;i&lt;st-&gt;subframeSize;i++)
            {
               spx_word16_t exci=exc[i];
               exc[i]= ADD16(ADD16(MULT16_16_Q15(QCONST16(.7f,15),exc[i]) , MULT16_16_Q15(QCONST16(.3f,15),st-&gt;voc_m1)),
                             SUB16(MULT16_16_Q15(Q15_ONE-MULT16_16_16(QCONST16(.85f,9),g),EXTRACT16(PSHR32(innov[i],SIG_SHIFT))),
                                   MULT16_16_Q15(MULT16_16_16(QCONST16(.15f,9),g),EXTRACT16(PSHR32(st-&gt;voc_m2,SIG_SHIFT)))
                                  ));
               st-&gt;voc_m1 = exci;
               st-&gt;voc_m2=innov[i];
               st-&gt;voc_mean = EXTRACT16(PSHR32(ADD32(MULT16_16(QCONST16(.8f,15),st-&gt;voc_mean), MULT16_16(QCONST16(.2f,15),exc[i])), 15));
               exc[i]-=st-&gt;voc_mean;
            }
         }

      }
   }
   
   ALLOC(interp_qlsp, st-&gt;lpcSize, spx_lsp_t);

   if (st-&gt;lpc_enh_enabled &amp;&amp; SUBMODE(comb_gain)&gt;0 &amp;&amp; !st-&gt;count_lost)
   {
      multicomb(st-&gt;exc-st-&gt;subframeSize, out, st-&gt;interp_qlpc, st-&gt;lpcSize, 2*st-&gt;subframeSize, best_pitch, 40, SUBMODE(comb_gain), stack);
      multicomb(st-&gt;exc+st-&gt;subframeSize, out+2*st-&gt;subframeSize, st-&gt;interp_qlpc, st-&gt;lpcSize, 2*st-&gt;subframeSize, best_pitch, 40, SUBMODE(comb_gain), stack);
   } else {
      SPEEX_COPY(out, &amp;st-&gt;exc[-st-&gt;subframeSize], st-&gt;frameSize);
   }
   
   /* If the last packet was lost, re-scale the excitation to obtain the same energy as encoded in ol_gain */
   if (st-&gt;count_lost) 
   {
      spx_word16_t exc_ener;
      spx_word32_t gain32;
      spx_word16_t gain;
      exc_ener = compute_rms16 (st-&gt;exc, st-&gt;frameSize);
      gain32 = PDIV32(ol_gain, ADD16(exc_ener,1));
#ifdef FIXED_POINT
      if (gain32 &gt; 32767)
         gain32 = 32767;
      gain = EXTRACT16(gain32);
#else
      if (gain32 &gt; 2)
         gain32=2;
      gain = gain32;
#endif
      for (i=0;i&lt;st-&gt;frameSize;i++)
      {
         st-&gt;exc[i] = MULT16_16_Q14(gain, st-&gt;exc[i]);
         out[i]=st-&gt;exc[i-st-&gt;subframeSize];
      }
   }

   /*Loop on subframes */
   for (sub=0;sub&lt;st-&gt;nbSubframes;sub++)
   {
      int offset;
      spx_word16_t *sp;
      spx_word16_t *exc;
      /* Offset relative to start of frame */
      offset = st-&gt;subframeSize*sub;
      /* Original signal */
      sp=out+offset;
      /* Excitation */
      exc=st-&gt;exc+offset;

      /* LSP interpolation (quantized and unquantized) */
      lsp_interpolate(st-&gt;old_qlsp, qlsp, interp_qlsp, st-&gt;lpcSize, sub, st-&gt;nbSubframes);

      /* Make sure the LSP&#39;s are stable */
      lsp_enforce_margin(interp_qlsp, st-&gt;lpcSize, LSP_MARGIN);

      /* Compute interpolated LPCs (unquantized) */
      lsp_to_lpc(interp_qlsp, ak, st-&gt;lpcSize, stack);

      /* Compute analysis filter at w=pi */
      {
         spx_word32_t pi_g=LPC_SCALING;
         for (i=0;i&lt;st-&gt;lpcSize;i+=2)
         {
            /*pi_g += -st-&gt;interp_qlpc[i] +  st-&gt;interp_qlpc[i+1];*/
            pi_g = ADD32(pi_g, SUB32(EXTEND32(ak[i+1]),EXTEND32(ak[i])));
         }
         st-&gt;pi_gain[sub] = pi_g;
      }
      
      iir_mem16(sp, st-&gt;interp_qlpc, sp, st-&gt;subframeSize, st-&gt;lpcSize, 
                st-&gt;mem_sp, stack);
      
      for (i=0;i&lt;st-&gt;lpcSize;i++)
         st-&gt;interp_qlpc[i] = ak[i];

   }

   if (st-&gt;highpass_enabled)
      highpass(out, out, st-&gt;frameSize, (st-&gt;isWideband?HIGHPASS_WIDEBAND:HIGHPASS_NARROWBAND)|HIGHPASS_OUTPUT, st-&gt;mem_hp);
   /*for (i=0;i&lt;st-&gt;frameSize;i++)
     printf (&quot;%d\n&quot;, (int)st-&gt;frame[i]);*/

   /* Tracking output level */
   st-&gt;level = 1+PSHR32(ol_gain,SIG_SHIFT);
   st-&gt;max_level = MAX16(MULT16_16_Q15(QCONST16(.99f,15), st-&gt;max_level), st-&gt;level);
   st-&gt;min_level = MIN16(ADD16(1,MULT16_16_Q14(QCONST16(1.01f,14), st-&gt;min_level)), st-&gt;level);
   if (st-&gt;max_level &lt; st-&gt;min_level+1)
      st-&gt;max_level = st-&gt;min_level+1;
   /*printf (&quot;%f %f %f %d\n&quot;, og, st-&gt;min_level, st-&gt;max_level, update);*/
   
   /* Store the LSPs for interpolation in the next frame */
   for (i=0;i&lt;st-&gt;lpcSize;i++)
      st-&gt;old_qlsp[i] = qlsp[i];

   /* The next frame will not be the first (Duh!) */
   st-&gt;first = 0;
   st-&gt;count_lost=0;
   st-&gt;last_pitch = best_pitch;
#ifdef FIXED_POINT
   st-&gt;last_pitch_gain = PSHR16(pitch_average,2);
#else
   st-&gt;last_pitch_gain = .25*pitch_average;   
#endif
   st-&gt;pitch_gain_buf[st-&gt;pitch_gain_buf_idx++] = st-&gt;last_pitch_gain;
   if (st-&gt;pitch_gain_buf_idx &gt; 2) /* rollover */
      st-&gt;pitch_gain_buf_idx = 0;

   st-&gt;last_ol_gain = ol_gain;

   return 0;
}

int nb_encoder_ctl(void *state, int request, void *ptr)
{
   EncState *st;
   st=(EncState*)state;     
   switch(request)
   {
   case SPEEX_GET_FRAME_SIZE:
      (*(spx_int32_t*)ptr) = st-&gt;frameSize;
      break;
   case SPEEX_SET_LOW_MODE:
   case SPEEX_SET_MODE:
      st-&gt;submodeSelect = st-&gt;submodeID = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_LOW_MODE:
   case SPEEX_GET_MODE:
      (*(spx_int32_t*)ptr) = st-&gt;submodeID;
      break;
#ifndef DISABLE_VBR
      case SPEEX_SET_VBR:
      st-&gt;vbr_enabled = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_VBR:
      (*(spx_int32_t*)ptr) = st-&gt;vbr_enabled;
      break;
   case SPEEX_SET_VAD:
      st-&gt;vad_enabled = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_VAD:
      (*(spx_int32_t*)ptr) = st-&gt;vad_enabled;
      break;
   case SPEEX_SET_DTX:
      st-&gt;dtx_enabled = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_DTX:
      (*(spx_int32_t*)ptr) = st-&gt;dtx_enabled;
      break;
   case SPEEX_SET_ABR:
      st-&gt;abr_enabled = (*(spx_int32_t*)ptr);
      st-&gt;vbr_enabled = st-&gt;abr_enabled!=0;
      if (st-&gt;vbr_enabled) 
      {
         spx_int32_t i=10;
         spx_int32_t rate, target;
         float vbr_qual;
         target = (*(spx_int32_t*)ptr);
         while (i&gt;=0)
         {
            speex_encoder_ctl(st, SPEEX_SET_QUALITY, &amp;i);
            speex_encoder_ctl(st, SPEEX_GET_BITRATE, &amp;rate);
            if (rate &lt;= target)
               break;
            i--;
         }
         vbr_qual=i;
         if (vbr_qual&lt;0)
            vbr_qual=0;
         speex_encoder_ctl(st, SPEEX_SET_VBR_QUALITY, &amp;vbr_qual);
         st-&gt;abr_count=0;
         st-&gt;abr_drift=0;
         st-&gt;abr_drift2=0;
      }
      
      break;
   case SPEEX_GET_ABR:
      (*(spx_int32_t*)ptr) = st-&gt;abr_enabled;
      break;
#endif /* #ifndef DISABLE_VBR */
#if !defined(DISABLE_VBR) &amp;&amp; !defined(DISABLE_FLOAT_API)
   case SPEEX_SET_VBR_QUALITY:
      st-&gt;vbr_quality = (*(float*)ptr);
      break;
   case SPEEX_GET_VBR_QUALITY:
      (*(float*)ptr) = st-&gt;vbr_quality;
      break;
#endif /* !defined(DISABLE_VBR) &amp;&amp; !defined(DISABLE_FLOAT_API) */
   case SPEEX_SET_QUALITY:
      {
         int quality = (*(spx_int32_t*)ptr);
         if (quality &lt; 0)
            quality = 0;
         if (quality &gt; 10)
            quality = 10;
         st-&gt;submodeSelect = st-&gt;submodeID = ((const SpeexNBMode*)(st-&gt;mode-&gt;mode))-&gt;quality_map[quality];
      }
      break;
   case SPEEX_SET_COMPLEXITY:
      st-&gt;complexity = (*(spx_int32_t*)ptr);
      if (st-&gt;complexity&lt;0)
         st-&gt;complexity=0;
      break;
   case SPEEX_GET_COMPLEXITY:
      (*(spx_int32_t*)ptr) = st-&gt;complexity;
      break;
   case SPEEX_SET_BITRATE:
      {
         spx_int32_t i=10;
         spx_int32_t rate, target;
         target = (*(spx_int32_t*)ptr);
         while (i&gt;=0)
         {
            speex_encoder_ctl(st, SPEEX_SET_QUALITY, &amp;i);
            speex_encoder_ctl(st, SPEEX_GET_BITRATE, &amp;rate);
            if (rate &lt;= target)
               break;
            i--;
         }
      }
      break;
   case SPEEX_GET_BITRATE:
      if (st-&gt;submodes[st-&gt;submodeID])
         (*(spx_int32_t*)ptr) = st-&gt;sampling_rate*SUBMODE(bits_per_frame)/st-&gt;frameSize;
      else
         (*(spx_int32_t*)ptr) = st-&gt;sampling_rate*(NB_SUBMODE_BITS+1)/st-&gt;frameSize;
      break;
   case SPEEX_SET_SAMPLING_RATE:
      st-&gt;sampling_rate = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_SAMPLING_RATE:
      (*(spx_int32_t*)ptr)=st-&gt;sampling_rate;
      break;
   case SPEEX_RESET_STATE:
      {
         int i;
         st-&gt;bounded_pitch = 1;
         st-&gt;first = 1;
         for (i=0;i&lt;st-&gt;lpcSize;i++)
            st-&gt;old_lsp[i]= DIV32(MULT16_16(QCONST16(3.1415927f, LSP_SHIFT), i+1), st-&gt;lpcSize+1);
         for (i=0;i&lt;st-&gt;lpcSize;i++)
            st-&gt;mem_sw[i]=st-&gt;mem_sw_whole[i]=st-&gt;mem_sp[i]=st-&gt;mem_exc[i]=0;
         for (i=0;i&lt;st-&gt;frameSize+st-&gt;max_pitch+1;i++)
            st-&gt;excBuf[i]=st-&gt;swBuf[i]=0;
         for (i=0;i&lt;st-&gt;windowSize-st-&gt;frameSize;i++)
            st-&gt;winBuf[i]=0;
      }
      break;
   case SPEEX_SET_SUBMODE_ENCODING:
      st-&gt;encode_submode = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_SUBMODE_ENCODING:
      (*(spx_int32_t*)ptr) = st-&gt;encode_submode;
      break;
   case SPEEX_GET_LOOKAHEAD:
      (*(spx_int32_t*)ptr)=(st-&gt;windowSize-st-&gt;frameSize);
      break;
   case SPEEX_SET_PLC_TUNING:
      st-&gt;plc_tuning = (*(spx_int32_t*)ptr);
      if (st-&gt;plc_tuning&gt;100)
         st-&gt;plc_tuning=100;
      break;
   case SPEEX_GET_PLC_TUNING:
      (*(spx_int32_t*)ptr)=(st-&gt;plc_tuning);
      break;
#ifndef DISABLE_VBR
   case SPEEX_SET_VBR_MAX_BITRATE:
      st-&gt;vbr_max = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_VBR_MAX_BITRATE:
      (*(spx_int32_t*)ptr) = st-&gt;vbr_max;
      break;
#endif /* #ifndef DISABLE_VBR */
   case SPEEX_SET_HIGHPASS:
      st-&gt;highpass_enabled = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_HIGHPASS:
      (*(spx_int32_t*)ptr) = st-&gt;highpass_enabled;
      break;

   /* This is all internal stuff past this point */
   case SPEEX_GET_PI_GAIN:
      {
         int i;
         spx_word32_t *g = (spx_word32_t*)ptr;
         for (i=0;i&lt;st-&gt;nbSubframes;i++)
            g[i]=st-&gt;pi_gain[i];
      }
      break;
   case SPEEX_GET_EXC:
      {
         int i;
         for (i=0;i&lt;st-&gt;nbSubframes;i++)
            ((spx_word16_t*)ptr)[i] = compute_rms16(st-&gt;exc+i*st-&gt;subframeSize, st-&gt;subframeSize);
      }
      break;
#ifndef DISABLE_VBR
   case SPEEX_GET_RELATIVE_QUALITY:
      (*(float*)ptr)=st-&gt;relative_quality;
      break;
#endif /* #ifndef DISABLE_VBR */
   case SPEEX_SET_INNOVATION_SAVE:
      st-&gt;innov_rms_save = (spx_word16_t*)ptr;
      break;
   case SPEEX_SET_WIDEBAND:
      st-&gt;isWideband = *((spx_int32_t*)ptr);
      break;
   case SPEEX_GET_STACK:
      *((char**)ptr) = st-&gt;stack;
      break;
   default:
      speex_warning_int(&quot;Unknown nb_ctl request: &quot;, request);
      return -1;
   }
   return 0;
}

int nb_decoder_ctl(void *state, int request, void *ptr)
{
   DecState *st;
   st=(DecState*)state;
   switch(request)
   {
   case SPEEX_SET_LOW_MODE:
   case SPEEX_SET_MODE:
      st-&gt;submodeID = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_LOW_MODE:
   case SPEEX_GET_MODE:
      (*(spx_int32_t*)ptr) = st-&gt;submodeID;
      break;
   case SPEEX_SET_ENH:
      st-&gt;lpc_enh_enabled = *((spx_int32_t*)ptr);
      break;
   case SPEEX_GET_ENH:
      *((spx_int32_t*)ptr) = st-&gt;lpc_enh_enabled;
      break;
   case SPEEX_GET_FRAME_SIZE:
      (*(spx_int32_t*)ptr) = st-&gt;frameSize;
      break;
   case SPEEX_GET_BITRATE:
      if (st-&gt;submodes[st-&gt;submodeID])
         (*(spx_int32_t*)ptr) = st-&gt;sampling_rate*SUBMODE(bits_per_frame)/st-&gt;frameSize;
      else
         (*(spx_int32_t*)ptr) = st-&gt;sampling_rate*(NB_SUBMODE_BITS+1)/st-&gt;frameSize;
      break;
   case SPEEX_SET_SAMPLING_RATE:
      st-&gt;sampling_rate = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_SAMPLING_RATE:
      (*(spx_int32_t*)ptr)=st-&gt;sampling_rate;
      break;
   case SPEEX_SET_HANDLER:
      {
         SpeexCallback *c = (SpeexCallback*)ptr;
         st-&gt;speex_callbacks[c-&gt;callback_id].func=c-&gt;func;
         st-&gt;speex_callbacks[c-&gt;callback_id].data=c-&gt;data;
         st-&gt;speex_callbacks[c-&gt;callback_id].callback_id=c-&gt;callback_id;
      }
      break;
   case SPEEX_SET_USER_HANDLER:
      {
         SpeexCallback *c = (SpeexCallback*)ptr;
         st-&gt;user_callback.func=c-&gt;func;
         st-&gt;user_callback.data=c-&gt;data;
         st-&gt;user_callback.callback_id=c-&gt;callback_id;
      }
      break;
   case SPEEX_RESET_STATE:
      {
         int i;
         for (i=0;i&lt;st-&gt;lpcSize;i++)
            st-&gt;mem_sp[i]=0;
         for (i=0;i&lt;st-&gt;frameSize + st-&gt;max_pitch + 1;i++)
            st-&gt;excBuf[i]=0;
      }
      break;
   case SPEEX_SET_SUBMODE_ENCODING:
      st-&gt;encode_submode = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_SUBMODE_ENCODING:
      (*(spx_int32_t*)ptr) = st-&gt;encode_submode;
      break;
   case SPEEX_GET_LOOKAHEAD:
      (*(spx_int32_t*)ptr)=st-&gt;subframeSize;
      break;
   case SPEEX_SET_HIGHPASS:
      st-&gt;highpass_enabled = (*(spx_int32_t*)ptr);
      break;
   case SPEEX_GET_HIGHPASS:
      (*(spx_int32_t*)ptr) = st-&gt;highpass_enabled;
      break;
      /* FIXME: Convert to fixed-point and re-enable even when float API is disabled */
#ifndef DISABLE_FLOAT_API
   case SPEEX_GET_ACTIVITY:
   {
      float ret;
      ret = log(st-&gt;level/st-&gt;min_level)/log(st-&gt;max_level/st-&gt;min_level);
      if (ret&gt;1)
         ret = 1;
      /* Done in a strange way to catch NaNs as well */
      if (!(ret &gt; 0))
         ret = 0;
      /*printf (&quot;%f %f %f %f\n&quot;, st-&gt;level, st-&gt;min_level, st-&gt;max_level, ret);*/
      (*(spx_int32_t*)ptr) = (int)(100*ret);
   }
   break;
#endif
   case SPEEX_GET_PI_GAIN:
      {
         int i;
         spx_word32_t *g = (spx_word32_t*)ptr;
         for (i=0;i&lt;st-&gt;nbSubframes;i++)
            g[i]=st-&gt;pi_gain[i];
      }
      break;
   case SPEEX_GET_EXC:
      {
         int i;
         for (i=0;i&lt;st-&gt;nbSubframes;i++)
            ((spx_word16_t*)ptr)[i] = compute_rms16(st-&gt;exc+i*st-&gt;subframeSize, st-&gt;subframeSize);
      }
      break;
   case SPEEX_GET_DTX_STATUS:
      *((spx_int32_t*)ptr) = st-&gt;dtx_enabled;
      break;
   case SPEEX_SET_INNOVATION_SAVE:
      st-&gt;innov_save = (spx_word16_t*)ptr;
      break;
   case SPEEX_SET_WIDEBAND:
      st-&gt;isWideband = *((spx_int32_t*)ptr);
      break;
   case SPEEX_GET_STACK:
      *((char**)ptr) = st-&gt;stack;
      break;
   default:
      speex_warning_int(&quot;Unknown nb_ctl request: &quot;, request);
      return -1;
   }
   return 0;
}
