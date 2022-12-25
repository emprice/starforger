import Two from 'two.js';
import { Path } from 'two.js/src/path';
import { Circle } from 'two.js/src/shapes/circle';

import dat from 'dat.gui';

import mandelagol from 'autogen/mandelagol.js';
import 'autogen/mandelagol.wasm';

interface IPointerArray {
  array: Float64Array;
  ptr: number;
}

interface ITransitWidgetConfig {
  p: number;
  gamma1: number;
  gamma2: number;
  playing: boolean;
}

export class TransitWidget {

  constructor(elem: HTMLElement) {

    const viewport = document.createElement('div');
    viewport.classList.add('transitWidgetViewport');

    const controls = document.createElement('div');
    controls.classList.add('transitWidgetControls');

    elem.appendChild(viewport);
    elem.appendChild(controls);

    this.gui = new dat.GUI({
      name: 'Transit parameters',
      hideable: false,
      closed: false,
      closeOnTop: false
    });
    this.gui.domElement = controls;

    this.config = {
      p: 0.1,
      gamma1: 0.05,
      gamma2: 0.01,
      playing: true
    };

    var pControl = this.gui.add(this.config, 'p', 0.01, 1, 0.01);
    var gamma1Control = this.gui.add(this.config, 'gamma1', 0, 1, 0.01);
    var gamma2Control = this.gui.add(this.config, 'gamma2', 0, 1, 0.01);
    var playControl = this.gui.add(this.config, 'playing');

    pControl.onChange((val: number) => {

      this.two.pause();
      this.two.frameCount = 0;
    });

    gamma1Control.onChange((val: number) => {

      this.two.pause();
      this.two.frameCount = 0;
    });

    gamma2Control.onChange((val: number) => {

      this.two.pause();
      this.two.frameCount = 0;
    });

    pControl.onFinishChange((val: number) => {

      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    gamma1Control.onFinishChange((val: number) => {

      const otherVal = gamma2Control.getValue();
      if (val + otherVal >= 1) {
        gamma1Control.setValue(1 - otherVal);
      }

      this.updateLimbDarkening();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    gamma2Control.onFinishChange((val: number) => {

      const otherVal = gamma1Control.getValue();
      if (val + otherVal >= 1) {
        gamma2Control.setValue(1 - otherVal);
      }

      this.updateLimbDarkening();
      this.updateTransit();
      if (this.config.playing) this.two.play();
    });

    playControl.onFinishChange((val: boolean) => {

      if (val) this.two.play();
      else this.two.pause();
    });

    const width = 800, height = 450;

    this.two = new Two({
      fullscreen: false,
      fitted: false,
      autostart: false,
      width: width,
      height: height,
      type: Two.Types.svg
    }).appendTo(viewport);

    this.two.renderer.domElement.setAttribute('preserveAspectRatio', 'xMidYMid meet');
    this.two.renderer.domElement.setAttribute('viewBox', `0 0 ${width} ${height}`);
    this.two.renderer.domElement.style.width = '100%';
    this.two.renderer.domElement.style.height = '100%';

    this.radius = this.two.width * 0.25;
    this.x = this.two.width * 0.5;
    this.y = this.two.height * 0.5;

    this.circle = this.two.makeCircle(this.x, this.y, this.radius);
    this.circle.stroke = 'black';
    this.circle.linewidth = 1;
    this.updateLimbDarkening();

    this.mandelagol = mandelagol();
    this.mandelagol.then((instance: any) => {

      this.z = TransitWidget.allocArray(instance, TransitWidget.nsamp);
      this.flux = TransitWidget.allocArray(instance, TransitWidget.nsamp);
    });

    this.curve = null;
    this.updateTransit();

    this.two.bind('update', (frame: number) => {
      const frac = (frame % 300) / 300;
      if (this.curve) {
        this.curve.ending = frac;
      }
    });
    this.two.play();
  }

  public finalize() {

    this.mandelagol.then((instance: any) => {

      TransitWidget.freeArray(instance, this.z);
      TransitWidget.freeArray(instance, this.flux);
    });
  }

  protected updateTransit() {

    this.mandelagol.then((instance: any) => {

      const zout = 1 + 2 * this.config.p;

      for (var i = 0; i < TransitWidget.nsamp; ++i) {

        this.z.array[i] = -zout + 2 * zout * (i / (TransitWidget.nsamp - 1));
      }

      instance._mandelagolArray(this.config.p, this.z.ptr,
        this.config.gamma1, this.config.gamma2, this.flux.ptr,
        TransitWidget.nsamp);

      var xy = new Array();

      for (var i = 0; i < TransitWidget.nsamp; ++i) {

        const z = this.z.array[i];
        const f = this.flux.array[i];

        const x = this.x + this.radius * z;
        const y = this.y + 10 * (1 - f) * this.radius;

        xy.push(x, y);
      }

      if (this.curve) this.curve.remove();
      this.curve = this.two.makePath(...xy);
      this.curve.noFill();
      this.curve.stroke = '#bf616a';
      this.curve.linewidth = 4;
      this.curve.curved = true;
      this.curve.closed = false;
    });
  }

  protected updateLimbDarkening() {

    var stops = new Array(), n = 25;

    for (var i = 0; i < n; ++i) {

      const mu1 = i / (n - 1);
      const intensity = 1. - this.config.gamma1 * mu1 -
        this.config.gamma2 * mu1 * mu1;
      const ix = (255 * intensity).toFixed(0);
      const color = `rgb(${ix}, ${ix}, ${ix})`;

      stops.push(new Two.Stop(mu1, color));
    }

    var gradient =
      this.two.makeRadialGradient(0, 0, this.radius, ...stops);
    gradient.units = 'userSpaceOnUse';
    this.circle.fill = gradient;
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

  static nsamp = 250;

  protected two: Two;
  protected gui: dat.GUI;

  protected x: number;
  protected y: number;
  protected radius: number;

  protected circle: Circle;
  protected curve: Path | null;

  protected config: ITransitWidgetConfig;

  protected z!: IPointerArray;
  protected flux!: IPointerArray;
  protected mandelagol: Promise<any>;
}

// vim: set ft=typescript:
