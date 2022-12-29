import Two from 'two.js';
import { Path } from 'two.js/src/path';
import { Circle } from 'two.js/src/shapes/circle';

import * as dat from 'dat.gui/src/dat/index.js';
import * as introJs from 'intro.js/intro.js';

import compute from 'autogen/compute.js';
import 'autogen/compute.wasm';

/// Wrapper for an array backed by WASM-accessible memory
interface IPointerArray {
  array: Float64Array;  ///< Javascript array
  ptr: number;          ///< WASM memory pointer
}

/// Configuration parameters for the widget
interface ITransitWidgetConfig {
  p: number;            ///< Radius ratio Rp/Rs
  gamma1: number;       ///< Linear limb darkening coefficient
  gamma2: number;       ///< Quadratic limb darkening coefficient
  Teff: number;         ///< Stellar effective temperature
  mass: number;         ///< Mass ratio Ms/Mp
  a: number;            ///< Semimajor axis a/Rs
  ecc: number;          ///< Orbital eccentricity
  omega: number;        ///< Argument of periastron
  playing: boolean;     ///< Is the simulation currently animating?
  explain: () => void;  ///< Run the intro tour again
}

export class TransitWidget {

  constructor(elem: HTMLElement) {

    // create the WASM instance
    this.compute = compute();

    // when the library is ready...
    this.compute.then((instance: any) => {

      // ... allocate all WASM-accessible arrays
      this.rgb   = TransitWidget.allocArray(instance, 3);
      this.t     = TransitWidget.allocArray(instance, TransitWidget.nsamp);
      this.r     = TransitWidget.allocArray(instance, TransitWidget.nsamp);
      this.theta = TransitWidget.allocArray(instance, TransitWidget.nsamp);
      this.x     = TransitWidget.allocArray(instance, TransitWidget.nsamp);
      this.y     = TransitWidget.allocArray(instance, TransitWidget.nsamp);
      this.flux  = TransitWidget.allocArray(instance, TransitWidget.nsamp);
    });

    // create the main widget container
    const viewport = document.createElement('div');
    viewport.classList.add('transitWidgetViewport');
    elem.appendChild(viewport);

    // add parameter controls
    this.gui = new dat.GUI({
      name: 'Parameters',
      hideable: false,
      closed: false,
      closeOnTop: false
    });

    // set up initial/default parameters
    this.config = {
      p: 0.2,
      gamma1: 0.5,
      gamma2: 0.01,
      Teff: 4000,
      mass: 100,
      a: 1000,
      ecc: 0.1,
      omega: 0.5 * Math.PI,
      playing: true,
      explain: () => {
        this.tour.setOption('dontShowAgain', false);
        this.two.pause();
        this.tour.start();
      }
    };

    // set up parameter groups
    var folder1 = this.gui.addFolder('Transit parameters');
    var folder2 = this.gui.addFolder('Orbit parameters');

    // set up transit parameter sliders
    var pControl      = folder1.add(this.config, 'p', 0.01, 1, 0.01);
    var gamma1Control = folder1.add(this.config, 'gamma1', 0, 1, 0.01);
    var gamma2Control = folder1.add(this.config, 'gamma2', 0, 1, 0.01);
    var teffControl   = folder1.add(this.config, 'Teff', 1000, 10000, 100);

    // set up orbit parameter sliders
    var massControl  = folder2.add(this.config, 'mass', 10, 1000, 10);
    var aControl     = folder2.add(this.config, 'a', 100, 10000, 100);
    var eccControl   = folder2.add(this.config, 'ecc', 0, 1, 0.01);
    var omegaControl = folder2.add(this.config, 'omega', 0, 2 * Math.PI, 0.05 * Math.PI);

    // set up global controls
    var playControl = this.gui.add(this.config, 'playing');
    var explainControl = this.gui.add(this.config, 'explain');

    var allControls = [pControl, gamma1Control, gamma2Control, teffControl,
      massControl, aControl, eccControl, omegaControl];

    for (var control of allControls) {

      control.onChange((val: number) => {

        this.two.pause();
        this.two.frameCount = 0;
      });
    }

    pControl.onFinishChange((val: number) => {

      this.updateOrbit();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    gamma1Control.onFinishChange((val: number) => {

      // clip to allowed value range gamma1 + gamma2 < 1
      const otherVal = gamma2Control.getValue();
      if (val + otherVal >= 1) {
        gamma1Control.setValue(1 - otherVal);
      }

      this.updateLimbDarkening();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    gamma2Control.onFinishChange((val: number) => {

      // clip to allowed value range gamma1 + gamma2 < 1
      const otherVal = gamma1Control.getValue();
      if (val + otherVal >= 1) {
        gamma2Control.setValue(1 - otherVal);
      }

      this.updateLimbDarkening();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    teffControl.onFinishChange((val: number) => {

      this.updateLimbDarkening();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    massControl.onFinishChange((val: number) => {

      this.updateOrbit();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    aControl.onFinishChange((val: number) => {

      this.updateOrbit();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    eccControl.onFinishChange((val: number) => {

      this.updateOrbit();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    omegaControl.onFinishChange((val: number) => {

      this.updateOrbit();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    playControl.onFinishChange((val: boolean) => {

      if (val) this.two.play();
      else this.two.pause();
    });

    // allow spacebar to toggle the play control
    document.addEventListener('keydown', (event) => {

      if (event.code === 'Space') {
        playControl.setValue(!this.config.playing);
      }
    });

    // widget size, arbitrarily chosen to be 16:9
    const width = 800, height = 450;

    this.two = new Two({
      fullscreen: false,
      fitted: false,
      autostart: false,
      width: width,
      height: height,
      type: Two.Types.svg
    }).appendTo(viewport);

    // grow to fill available space while keeping aspect ratio; two.js doesn't
    // set these, but we can override them manually
    this.two.renderer.domElement.setAttribute('preserveAspectRatio', 'xMidYMid meet');
    this.two.renderer.domElement.setAttribute('viewBox', `0 0 ${width} ${height}`);
    this.two.renderer.domElement.style.width = '100%';
    this.two.renderer.domElement.style.height = '100%';

    // set up object sizes and center canvas position
    this.radius = this.two.width * 0.25;
    this.x0 = this.two.width * 0.5;
    this.y0 = this.two.height * 0.5;

    // create the stellar disk
    this.disk = this.two.makeCircle(this.x0, this.y0, this.radius);
    this.disk.noStroke();
    this.disk.noFill();

    // add the limb darkening
    this.compColor = '';
    this.updateLimbDarkening();

    // compute the initial orbit
    this.updateOrbit();

    // initialize the transit curve
    this.transit = null;
    this.updateTransit();

    this.two.bind('update', (frame: number) => {

      // frame loop draws transit moving across the stellar disk
      const idx = frame % TransitWidget.nsamp;
      if (this.transit) {
        this.transit.ending = this.t.array[idx] /
          this.t.array[TransitWidget.nsamp-1];
      }
    });

    // set up our nifty how-to tour
    this.tour = introJs().setOptions({
      disableInteraction: true,
      showStepNumbers: true,
      dontShowAgain: true,
      steps: [
        {
          title: 'Hello, world!',
          intro: 'Welcome to <b>STARFORGER</b>, the fanciest real-time transit light curve simulator on the Web. Possibly the Universe.'
        },
        {
          element: this.gui.domElement,
          intro: 'You can use these controls to change parameters that affect the planet\'s transit.',
          position: 'left'
        },
        {
          element: folder1.domElement,
          intro: 'These controls change the relative planet size, limb darkening parameters, and stellar effective temperature.',
          position: 'left'
        },
        {
          element: folder2.domElement,
          intro: 'These controls change the physical parameters of the star-planet system, thereby modifying the orbit.',
          position: 'left'
        },
        {
          title: 'All done!',
          intro: 'To pause the simulation, press the spacebar at any time.'
        }
      ]
    });

    // open folders right before they are introduced
    this.tour.onbeforechange((elem: any) => {
      if (elem === folder1.domElement) {
        folder1.open();
      } else if (elem === folder2.domElement) {
        folder2.open();
      }
    });

    // start the tour as soon as possible; when everything is ready,
    // start the animation
    this.compute.then((instance) => {

      this.two.update();
      this.tour.start();
      this.tour.oncomplete(() => {
        this.two.play();
      });
    });
  }

  public finalize() {

    // when the library is ready...
    this.compute.then((instance: any) => {

      // ... free the allocated arrays
      TransitWidget.freeArray(instance, this.rgb);
      TransitWidget.freeArray(instance, this.t);
      TransitWidget.freeArray(instance, this.r);
      TransitWidget.freeArray(instance, this.theta);
      TransitWidget.freeArray(instance, this.x);
      TransitWidget.freeArray(instance, this.y);
      TransitWidget.freeArray(instance, this.flux);
    });
  }

  protected updateOrbit() {

    // when the library is ready...
    this.compute.then((instance: any) => {

      const m = 1, M = this.config.mass;
      const zout = 1 + 2 * this.config.p;

      // ... re-calculate the orbit times and positions
      instance._orbitArray(M, m, this.config.a, this.config.ecc,
        this.config.omega, -zout, zout, this.t.ptr, this.r.ptr,
        this.theta.ptr, this.x.ptr, this.y.ptr, TransitWidget.nsamp);
    });
  }

  protected updateTransit() {

    // when the library is ready...
    this.compute.then((instance: any) => {

      // ... re-calculate the transit fluxes
      instance._mandelagolArray(this.config.p, this.x.ptr, this.config.gamma1,
        this.config.gamma2, this.flux.ptr, TransitWidget.nsamp);

      var xy = new Array();

      for (var i = 0; i < TransitWidget.nsamp; ++i) {

        const z = this.x.array[i];
        const f = this.flux.array[i];

        // screen position from (z, flux)
        // NOTE: the 10 here is not magic, it makes the transit deeper for
        // plotting on the screen
        const x = this.x0 + this.radius * z;
        const y = this.y0 + 10 * (1 - f) * this.radius;

        xy.push(x, y);
      }

      // re-draw the curve
      if (this.transit) this.transit.remove();
      this.transit = this.two.makePath(...xy);
      this.transit.noFill();
      this.transit.stroke = this.compColor;
      this.transit.linewidth = 4;
      this.transit.curved = true;
      this.transit.closed = false;
    });
  }

  protected updateLimbDarkening() {

    const clip = (val: number) => {
      return Math.max(Math.min(val, 1), 0);
    };

    // when the library is ready...
    this.compute.then((instance) => {

      // ... compute the blackbody color values for a given Teff
      instance._blackbodyToRgb(this.config.Teff, this.rgb.ptr);

      var stops = new Array(), n = 25;

      for (var i = 0; i < n; ++i) {

        const mu1 = i / (n - 1);
        const intensity = 1. - this.config.gamma1 * mu1 -
          this.config.gamma2 * mu1 * mu1;

        const ir = (255 * clip(intensity * this.rgb.array[0])).toFixed(0);
        const ig = (255 * clip(intensity * this.rgb.array[1])).toFixed(0);
        const ib = (255 * clip(intensity * this.rgb.array[2])).toFixed(0);

        // add a gradient stop for this color
        const color = `rgb(${ir}, ${ig}, ${ib})`;
        stops.push(new Two.Stop(mu1, color));
      }

      // add the gradient to the stellar disk fill
      var gradient =
        this.two.makeRadialGradient(0, 0, this.radius, ...stops);
      gradient.units = 'userSpaceOnUse';
      this.disk.fill = gradient;

      let hue: number;
      const ir = clip(this.rgb.array[0]);
      const ig = clip(this.rgb.array[1]);
      const ib = clip(this.rgb.array[2]);

      // calculate the hue for the non-darkened stellar color
      // https://stackoverflow.com/a/23094494/1552418
      if (ir > Math.max(ig, ib)) {
        hue = 60 * (ig - ib) / (ir - Math.min(ig, ib));
      } else if (ig > Math.max(ir, ib)) {
        hue = 60 * (2 + (ib - ir) / (ig - Math.min(ir, ib)));
      } else {
        hue = 60 * (4 + (ir - ig) / (ib - Math.min(ir, ig)));
      }

      // make a nice complementary color by rotating 180 degrees
      hue = (((hue + 180) % 360) + 360) % 360;
      const h = hue.toFixed(0);
      this.compColor = `hsl(${h}deg 50% 50%)`;
    });
  }

  static allocArray(instance: any, n: number): IPointerArray {

    // https://stackoverflow.com/a/49448982/1552418
    const nb = Float64Array.BYTES_PER_ELEMENT;
    var offset = instance._malloc(n * nb);
    instance.HEAPF64.set(new Float64Array(n), offset / nb);

    return {
      array: instance.HEAPF64.subarray(offset / nb, offset / nb + n),
      ptr: offset
    };
  }

  static freeArray(instance: any, array: IPointerArray) {

    instance._free(array.ptr);
    array.ptr = 0;
  }

  static nsamp = 300;   ///< Number of transit samples, rendered 1/frame

  protected two: Two;
  protected gui: typeof dat.GUI;
  protected tour: typeof introJs;

  protected x0: number;     ///< Horizontal center of canvas
  protected y0: number;     ///< Vertical center of canvas
  protected radius: number; ///< Stellar disk size

  protected disk: Circle;           ///< Stellar disk object
  protected transit: Path | null;   ///< Transit light curve object
  protected compColor: string;      ///< Color for the transit stroke

  protected config: ITransitWidgetConfig;

  protected rgb!: IPointerArray;    ///< Array of single (r, g, b) triple

  protected t!: IPointerArray;      ///< Times of transit and orbit samples
  protected r!: IPointerArray;      ///< Orbital distance from planet to star
  protected theta!: IPointerArray;  ///< Orbital angle around ellipse

  protected x!: IPointerArray;      ///< Projection of orbit along x
  protected y!: IPointerArray;      ///< Projection of orbit along y

  protected flux!: IPointerArray;   ///< Transit flux

  /** Compute kernel object wrapped by promise, so it cannot be used before
   *  it completely loads */
  protected compute: Promise<any>;
}

// vim: set ft=typescript:
