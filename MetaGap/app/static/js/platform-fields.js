(function () {
  "use strict";

  function toArray(nodeList) {
    return Array.prototype.slice.call(nodeList || []);
  }

  function findPlatformSelector() {
    return document.querySelector("[data-platform-selector]");
  }

  function findPlatformScopes() {
    return toArray(document.querySelectorAll("[data-platform-scope]"));
  }

  function toggleScope(wrapper, isActive) {
    if (!wrapper) {
      return;
    }

    const controls = toArray(
      wrapper.querySelectorAll("input, select, textarea, button")
    );

    if (isActive) {
      wrapper.classList.remove("d-none");
      controls.forEach(function (element) {
        element.disabled = false;
      });
    } else {
      wrapper.classList.add("d-none");
      controls.forEach(function (element) {
        element.disabled = true;
      });
    }
  }

  function applyPlatformSelection(selector, wrappers) {
    if (!selector) {
      return;
    }

    const selected = selector.value;

    wrappers.forEach(function (wrapper) {
      const scope = wrapper.getAttribute("data-platform-scope");
      const isActive = !scope || scope === selected;
      toggleScope(wrapper, isActive);
    });
  }

  function init() {
    const selector = findPlatformSelector();
    if (!selector) {
      return;
    }

    const wrappers = findPlatformScopes();
    const onChange = function () {
      applyPlatformSelection(selector, wrappers);
    };

    selector.addEventListener("change", onChange);
    selector.addEventListener("input", onChange);

    applyPlatformSelection(selector, wrappers);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
